import os
# from myontology import *

def similarity2csv(v1: list, v2: list=None, sim_type='similarity', file=None, cartesian=True):
    str2write = 'Entity,Entity,Similarity\r\n'
    if v2 is None:
        for i in v1:
            for j in v1:
                str2write += str(i) + ',' + str(j) + ','
                if sim_type == 'similarity':
                    str2write += str(i.similarity(j)) + '\r\n'
                elif sim_type == 'taxonomic_similarity':
                    str2write += str(i.taxonomic_similarity(j)) + '\r\n'
                elif sim_type == 'similarity_neighbors':
                    str2write += str(i.similarity_neighbors(j)) + '\r\n'
    else:
        if cartesian:
            for i in v1:
                for j in v2:
                    str2write += str(i) + ',' + str(j) + ','
                    if sim_type == 'similarity':
                        str2write += str(i.similarity(j)) + '\r\n'
                    elif sim_type == 'taxonomic_similarity':
                        str2write += str(i.taxonomic_similarity(j)) + '\r\n'
                    elif sim_type == 'similarity_neighbors':
                        str2write += str(i.similarity_neighbors(j)) + '\r\n'
        else:
            for i, j in zip(v1, v2):
                str2write += str(i) + ',' + str(j) + ','
                if sim_type == 'similarity':
                    str2write += str(i.similarity(j)) + '\r\n'
                elif sim_type == 'taxonomic_similarity':
                    str2write += str(i.taxonomic_similarity(j)) + '\r\n'
                elif sim_type == 'similarity_neighbors':
                    str2write += str(i.similarity_neighbors(j)) + '\r\n'
    if file is None:
        with open('results/' + sim_type + '.csv', 'wt') as f:
            f.write(str2write)
    else:
        with open(file, 'wt') as f:
            f.write(str2write)

# calculate similarity matrix
def similarity_calculation(v1: list, v2: list, sim_type='similarity'):
    if v2 is None:
        sim_mat = np.zeros((len(v1), len(v1)))
        for i, u in enumerate(v1):
            for j, v in enumerate(v1):
                if sim_type == 'similarity':
                    sim_mat[i, j] = u.similarity(v)
                elif sim_type == 'taxonomic_similarity':
                    sim_mat[i, j] = u.taxonomic_similarity(v)
                elif sim_type == 'similarity_neighbors':
                    sim_mat[i, j] = u.similarity_neighbors(v)
    else:
        for i, u in enumerate(v1):
            for j, v in enumerate(v2):
                if sim_type == 'similarity':
                    sim_mat[i, j] = u.similarity(v)
                elif sim_type == 'taxonomic_similarity':
                    sim_mat[i, j] = u.taxonomic_similarity(v)
                elif sim_type == 'similarity_neighbors':
                    sim_mat[i, j] = u.similarity_neighbors(v)
    
    return sim_mat

# save sim_mat as sim_mat.npy file in data directory
def save_sim_mat(sim_mat):
    np.save('data/sim_mat.npy', sim_mat)

# read sim_mat from sim_mat.npy file in data directory
def read_sim_mat():
    return np.load('data/sim_mat.npy')

# --------------- METIS ----------------

# write similarity matrix to a graph file of METIS's format
def sim_mat_to_graph(sim_mat, output_path='results/sim_mat.gp'):
    sim_mat = np.array(sim_mat)
    n = str(sim_mat.shape[0])  # the number of vertices
    fmt = '001'  # the format of graph
    ncon = '0'  # the number of multi-constraints of vertices

    m = 0  # the number of edges
    graph_structure = ''
    threshold = 0.3  # only consider edges with weights greater than threshold
    for i, row in enumerate(sim_mat):
        v_w_list = []
        for v, w in enumerate(row):
            if w >= threshold and i != v:
                v_w_list.append(str(v+1))
                v_w_list.append(str(int(w*1e6)))  # weights in graph file must be integer, so scale 10 times
                m += 1
        graph_structure += ' '.join(v_w_list) + '\n'
    m = int(m / 2)  # rule out duplicated edges
    string_to_write = n + ' ' + str(m) + ' ' + fmt + ' ' + ncon + '\n' + graph_structure

    with open(output_path, 'wt') as f:
        f.write(string_to_write)

# extract entities from entities csv file into groups according to METIS group file (each group is also a csv file containing entities)
def group_entities_into_csvfiles(partition_file, entities_to_group='results/entities_from_query.csv'):
    df = pd.read_csv(entities_to_group, sep=',', error_bad_lines=False, encoding='cp1252')
    entities = df['Unique_Entities'].to_list()
    groups = {}  # group dictionary, each key has all entities from one group
    with open(partition_file, 'rt') as f:
        lines = f.readlines()
        lines = [line[0] for line in lines]
        for i, line in enumerate(lines):
            if line not in groups.keys():
                groups[line] = []
            groups[line].append(entities[i])
    for k in groups.keys():
        utilities.write_entities_to_csv(groups[k], 'results/metis-results/group' + k + '.csv')

# calculate conductance of communities (inter-cluster conductance)
def conductance_comm(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        df = pd.read_csv(community_dir + '/' + cf, sep=',', error_bad_lines=False, encoding='cp1252')
        entity_set = set(df['Unique_Entities'].to_list())
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    conductance_clusters = []
    for comm in communities:
        con = nx.algorithms.cuts.conductance(G, comm)
        conductance_clusters.append(con)
    
    conductance = 1 - np.max(conductance_clusters)

    return conductance

# calculate coverage of communities
def coverage_comm(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        df = pd.read_csv(community_dir + '/' + cf, sep=',', error_bad_lines=False, encoding='cp1252')
        entity_set = set(df['Unique_Entities'].to_list())
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    # coverage
    coverage = nx.algorithms.community.quality.coverage(G, communities)
    
    return coverage

# calculate modularity of communities
def modularity_comm(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        df = pd.read_csv(community_dir + '/' + cf, sep=',', error_bad_lines=False, encoding='cp1252')
        entity_set = set(df['Unique_Entities'].to_list())
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    # coverage
    modularity = nx.community.modularity(G, communities)
    
    return modularity

# calculate performance of communities
def performance_comm(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        df = pd.read_csv(community_dir + '/' + cf, sep=',', error_bad_lines=False, encoding='cp1252')
        entity_set = set(df['Unique_Entities'].to_list())
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    # coverage
    performance = nx.algorithms.community.quality.performance(G, communities)
    
    return performance

# calculate normalized total cut of communities
def normalized_total_cut_comm(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        df = pd.read_csv(community_dir + '/' + cf, sep=',', error_bad_lines=False, encoding='cp1252')
        entity_set = set(df['Unique_Entities'].to_list())
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    total_cut = 0
    for i in range(len(communities)):
        for j in range(i + 1, len(communities)):
            total_cut += nx.cut_size(G, communities[i], communities[j], weight='weight')
    normalized_total_cut = total_cut / nx.volume(G, G, weight='weight')
    
    return normalized_total_cut


# --------------- semEP ----------------

# generate all graph files that are mandatory for semEP
def generate_bigraph(vertices1: list, vertices2: list, edges: list, outdir='/mnt/c/Users/SongZ/Downloads/repositories/semep/test/p4lucat/'):
    """
    Generate all graph files that are mandatory for semEP
    :param vertices1: a list contains elements of type MyOWLLogicalEntity
    :param vertices2: a list contains elements of type MyOWLLogicalEntity
    :param edges: a list conatins tuples, which describe the edge between vertices
                    for example: [(vertex in vertices1, vertex in vertices2, weight of edge)]
    :param outdir: the directory of output files
    
    :return: None
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # path of all necessary files
    vertices1_file = outdir + 'vertices1.txt'
    vertices1_simmat_file = outdir + 'vertices1_simmat.txt'
    vertices2_file = outdir + 'vertices2.txt'
    vertices2_simmat_file = outdir + 'vertices2_simmat.txt'
    bipartite_graph_file = outdir + 'bigraph.txt'

    with open(vertices1_file, 'wt') as v1f:
        v1f.write('{:d}\n'.format(len(vertices1)))
        for v1 in vertices1[:-1]:
            v1f.write('{}\n'.format(str(v1)))
        v1f.write('{}'.format(str(vertices1[-1])))
    
    with open(vertices1_simmat_file, 'wt') as v1sf:
        v1sf.write('{:d}\n'.format(len(vertices1)))
        for v11 in vertices1:
            for v12 in vertices1:
                sim = v11.similarity(v12)
                v1sf.write('{:.6f} '.format(sim))
                # v1sf.write('{:.6f} '.format(1.0))
                # v1sf.write('{:.6f} '.format(v11.taxonomic_similarity(v12)))
                print(str(v11), str(v12), sim)
            v1sf.write('\n')
    
    with open(vertices2_file, 'wt') as v2f:
        v2f.write('{:d}\n'.format(len(vertices2)))
        for v2 in vertices2[:-1]:
            v2f.write('{}\n'.format(str(v2)))
        v2f.write('{}'.format(str(vertices2[-1])))
    
    with open(vertices2_simmat_file, 'wt') as v2sf:
        v2sf.write('{:d}\n'.format(len(vertices2)))
        for v21 in vertices2:
            for v22 in vertices2:
                sim = v21.similarity(v22)
                v2sf.write('{:.6f} '.format(sim))
                # v2sf.write('{:.6f} '.format(1.0))
                # v2sf.write('{:.6f} '.format(v21.taxonomic_similarity(v22)))
                print(str(v21), str(v22), sim)
            v2sf.write('\n')
    
    with open(bipartite_graph_file, 'wt') as bgf:
        bgf.write('{:d}\n'.format(len(edges)))
        for v1, v2, w in edges:
            if w >= 0.0:
                bgf.write('{}\t{}\t{:.8f}\n'.format(str(v1), str(v2), w))

# calculate conductance of communities (inter-cluster conductance)
def comm_conductance(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        entity_set = set()
        with open(community_dir + '/' + cf) as cfile:
            lines = cfile.readlines()
            for line in lines:
                entity_set.update(line.split('\t')[:2])
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    conductance_clusters = []
    for comm in communities:
        con = nx.algorithms.cuts.conductance(G, comm)
        conductance_clusters.append(con)
    
    conductance = 1 - np.max(conductance_clusters)

    return conductance

# calculate coverage of communities
def comm_coverage(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        entity_set = set()
        with open(community_dir + '/' + cf) as cfile:
            lines = cfile.readlines()
            for line in lines:
                entity_set.update(line.split('\t')[:2])
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    # coverage
    coverage = nx.algorithms.community.quality.coverage(G, communities)
    print(coverage)

# calculate modularity of communities
def comm_modularity(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        entity_set = set()
        with open(community_dir + '/' + cf) as cfile:
            lines = cfile.readlines()
            for line in lines:
                entity_set.update(line.split('\t')[:2])
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    # coverage
    modularity = nx.community.modularity(G, communities)
    
    return modularity

# calculate performance of communities
def comm_performance(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        entity_set = set()
        with open(community_dir + '/' + cf) as cfile:
            lines = cfile.readlines()
            for line in lines:
                entity_set.update(line.split('\t')[:2])
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    # coverage
    performance = nx.algorithms.community.quality.performance(G, communities)
    
    return performance

# calculate normalized total cut of communities
def comm_normalized_total_cut(entities, sim_mat, community_dir, threshold=0.0):
    num_nodes = sim_mat.shape[0]

    # create graph from similarity matrix
    G = nx.Graph()
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if sim_mat[i, j] >= threshold:
                G.add_edge(i, j, weight=sim_mat[i, j])
    
    ### get communities from files ###
    communities = []

    com_files = os.listdir(community_dir)
    com_num = len(com_files)  # number of communities
    for i, cf in enumerate(com_files):
        entity_set = set()
        with open(community_dir + '/' + cf) as cfile:
            lines = cfile.readlines()
            for line in lines:
                entity_set.update(line.split('\t')[:2])
        print(f'Community {i}')
        print('---------------------------------------------')
        comm = [entities.index(e) for e in entity_set]
        communities.append(comm)
    
    total_cut = 0
    for i in range(len(communities)):
        for j in range(i + 1, len(communities)):
            total_cut += nx.cut_size(G, communities[i], communities[j], weight='weight')
    normalized_total_cut = total_cut / nx.volume(G, G, weight='weight')
    
    return normalized_total_cut

# call it for plotting radar chart
def radar_factory(num_vars, frame='circle'):
    """
    Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = 'radar'
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta

