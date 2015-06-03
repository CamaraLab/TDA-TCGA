"""
Copyright 2015, Pablo G. Camara, Columbia University
"""

import networkx
import multiprocessing
import json
import numpy
import numpy.linalg
import random
import requests
import scipy.stats
import scipy.cluster.hierarchy as sch
import pickle
import pylab
import copy_reg
import types


"""
Taken from https://gist.github.com/fiatmoney/1086393
"""


def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    if func_name.startswith('__') and not func_name.endswith('__'):
        cls_name = cls.__name__.lstrip('_')
        func_name = '_' + cls_name + func_name
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.__mro__:
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


"""
GLOBAL METHODS
"""


def ParseAyasdiGraph(lab, source, user, password, name):
    """
    Parses Ayasdi graph given by the source ID and lab ID, and stores as name.gexf and name.pickle. user and
    password specify Ayasdi login credentials.
    """
    headers = {"Content-type": "application/json"}
    session = requests.Session()
    session.post('https://core.ayasdi.com/login', data={'username': user, 'passphrase': password})
    r = session.get('https://core.ayasdi.com/v0/sources/' + source + '/networks/' + lab)
    sp = json.loads(r.content)
    rows = [int(x['id']) for x in sp['nodes']]
    dic2 = {}
    for i in rows:
        payload = {"network_nodes_descriptions": [{"network_id": lab, "node_ids": [i]}]}
        r = session.post('https://core.ayasdi.com/v0/sources/' + source + '/retrieve_row_indices',
                         data=json.dumps(payload), headers=headers)
        dic2[i] = json.loads(r.content)['row_indices']
    with open(name + '.pickle', 'wb') as handle3:
        pickle.dump(dic2, handle3)
    rowcount = []
    with open(name + '.gexf', 'w') as g:
        g.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        g.write('<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">\n')
        g.write('\t<graph mode="static" defaultedgetype="undirected">\n')
        g.write('\t\t<nodes>\n')
        for nod in sp['nodes']:
            g.write('\t\t\t<node id="' + str(nod['id']) + '" label="' + str(nod['row_count']) + '" />\n')
            rowcount.append(float(nod['row_count']))
        g.write('\t\t</nodes>\n')
        g.write('\t\t<edges>\n')
        for n5, edg in enumerate(sp['links']):
            g.write('\t\t\t<edge id="' + str(n5) + '" source="' + str(edg['from']) + '" target="' + str(edg['to'])
                    + '" />\n')
        g.write('\t\t</edges>\n')
        g.write('\t</graph>\n')
        g.write('</gexf>\n')


def benjamini_hochberg(pvalues):
    """
    Benjamini-Hochberg adjusted p-values for multiple testing. Consistent with R
    """
    pvalues = numpy.array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = numpy.empty(n)
    values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in xrange(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return list(new_pvalues)


def is_number(s):
    """
    Checks whether a string can be converted into a float
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def hierarchical_clustering(mat, method='median'):
    """
    Performs hierarchical clustering based on distance matrix mat and method
    """
    D = numpy.array(mat)
    fig = pylab.figure(figsize=(8, 8))
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])
    Y = sch.linkage(D, method=method)
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2 = fig.add_axes([0.3, 0.71, 0.6, 0.2])
    Y = sch.linkage(D, method=method)
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1, :]
    D = D[:, idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.get_cmap('jet_r'))
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.6])
    pylab.colorbar(im, cax=axcolor)
    pylab.show()
    return Z1


def _flat_connectivity(dicgen, pl, adj, ind=1):
    """
    Computes order ind connectivity of column genis
    """
    ex = []
    for uu in pl:
        ex.append(dicgen[uu])
    ex = numpy.array(ex)
    cor = float(len(pl))/float(len(pl)-1)
    return cor*numpy.dot(numpy.transpose(ex),
                         numpy.dot(numpy.linalg.matrix_power(adj, ind), ex))


"""
CLASSES
"""


class UnrootedGraph(object):
    """
    Graph: class initialized with the common name of a .gexf and .pickle file, and the name of a tab separated file
     containing the raw data, where the first row and first column contains ID's
    """
    def __init__(self, name, table, shift=None, log2=False, posgl=False):
        """
        Initialization function. name is the common name of .gexf and .pickle files. table is the name of the file
        containing the raw data. It ignores the first shift columns.
        """
        self.name = name
        self.g = networkx.read_gexf(name + '.gexf')
        self.gl = list(networkx.connected_component_subgraphs(self.g))[0]
        self.pl = self.gl.nodes()
        self.adj = numpy.array(networkx.to_numpy_matrix(self.gl, nodelist=self.pl))
        self.log2 = log2
        if not posgl:
            try:
                self.posgl = networkx.graphviz_layout(self.gl, 'sfdp', '-Goverlap=false -GK=0.1')
                self.posg = networkx.graphviz_layout(self.g, 'sfdp', '-Goverlap=false -GK=0.1')
            except:
                self.posgl = networkx.spring_layout(self.gl)
                self.posg = networkx.spring_layout(self.g)
            with open(name + '.posgl', 'w') as handler:
                pickle.dump(self.posgl, handler)
            with open(name + '.posg', 'w') as handler:
                pickle.dump(self.posg, handler)
        else:
            with open(name + '.posgl', 'r') as handler:
                self.posgl = pickle.load(handler)
            with open(name + '.posg', 'r') as handler:
                self.posg = pickle.load(handler)
        with open(name + '.pickle', 'r') as handler:
            self.dic = pickle.load(handler)
        with open(table, 'r') as f:
            self.dicgenes = {}
            self.geneindex = {}
            for n, line in enumerate(f):
                sp = numpy.array(line[:-1].split('\t'))
                if shift is None:
                    if n == 0:
                        poi = sp
                    else:
                        poi = []
                        for mji in sp:
                            if is_number(mji):
                                poi.append(mji)
                            else:
                                poi.append(0.0)
                elif type(shift) == int:
                    poi = sp[shift:]
                else:
                    poi = sp[shift]
                if n == 0:
                    for u, q in enumerate(poi):
                        self.dicgenes[q] = []
                        self.geneindex[u] = q
                else:
                    for u, q in enumerate(poi):
                        self.dicgenes[self.geneindex[u]].append(float(q))

    def get_gene(self, genin, permut=False, ignore_log=False):
        """
        Returns a dictionary that asigns to each node ID the average value of the column genin in the raw table, for
         the rows contained in that column. If permut = True performs a reshuffling of rows. If log2 = True assumes
          values in the format log2(1+x) when computing averages.
        """
        if type(genin) != list:
            genin = [genin]
        genecolor = {}
        for i in self.dic.keys():
            genecolor[str(i)] = 0.0
        for mju in genin:
            if mju is None:
                for i in sorted(self.dic.keys()):
                    genecolor[str(i)] = 0.0
            else:
                geys = self.dicgenes[mju]
                kis = range(len(geys))
                if permut:
                    random.shuffle(kis)
                for i in sorted(self.dic.keys()):
                    pol = 0.0
                    if self.log2 and not ignore_log:
                        for j in self.dic[i]:
                                pol += (numpy.power(2, float(geys[kis[j]]))-1.0)
                        pol = numpy.log2(1.0+(pol/float(len(self.dic[i]))))
                    else:
                        for j in self.dic[i]:
                                pol += float(geys[kis[j]])
                        pol /= float(len(self.dic[i]))
                    genecolor[str(i)] += pol
        tol = sum(genecolor.values())
        if tol > 0.0:
            for ll in genecolor.keys():
                genecolor[ll] = genecolor[ll]/tol
        return genecolor, tol

    def connectivity(self, genis, permut=False, ind=1):
        """
        Computes order ind connectivity of column genis
        """
        dicgen = self.get_gene(genis, permut)[0]
        return _flat_connectivity(dicgen, self.pl, self.adj, ind)

    def connectivity_pvalue(self, genis, ind=1, n=500, c=2):
        """
        Computes p-value of order ind connectivity of column genis performing permutation test with n permutations
        """
        if c > 1:
            pool = multiprocessing.Pool(processes=c, maxtasksperchild=1000)
            conn_async = [pool.apply_async(_flat_connectivity,
                                           [self.get_gene(genis, True)[0], self.pl, self.adj, ind]) for _ in range(n)]
            conn = [r.get() for r in conn_async]
            pool.close()
            pool.join()
        else:
            conn = [self.connectivity(genis, True, ind) for _ in range(n)]
        return float(sum(i > self.connectivity(genis) for i in conn))/float(n)

    def delta(self, genis):
        """
        Computes coefficient of variation of column genis in the graph
        """
        dicgen = self.get_gene(genis)[0]
        mi = [dicgen[node] for node in self.pl]
        if numpy.mean(mi) > 0.0:
            return numpy.std(mi)/numpy.mean(mi)
        else:
            return 0.0

    def expr(self, genis):
        """
        Computes fraction of expressed nodes for column genis
        """
        dicgen = self.get_gene(genis)[0]
        mi = [dicgen[node] for node in self.pl]
        return float(len(mi) - mi.count(0.0))/float(len(mi))

    def save(self, n=500, c=2):
        """
        Computes and stores all statistics in file name
        """
        pol = []
        with open(self.name + '.genes.txt', 'w') as ggg:
            ggg.write('Gene\tNodes\tDelta\tConnectivity\tp_value\tBH p-value\n')
            for gi in sorted(self.dicgenes.keys()):
                po = self.expr(gi)
                if po > 0:
                    pol.append(self.connectivity_pvalue(gi, ind=1, n=n, c=c))
            por = benjamini_hochberg(pol)
            mj = 0
            for gi in sorted(self.dicgenes.keys()):
                po = self.expr(gi)
                if po > 0:
                    ggg.write(gi + '\t' + str(po) + '\t' + str(self.delta(gi)) + '\t' +
                              str(self.connectivity(gi, ind=1)) + '\t' +
                              str(pol[mj]) + '\t' + str(por[mj]) + '\n')
                    mj += 1

    def JSD_matrix(self, lista):
        """
        Computes the Jensen-Shannon distance matrix of the list of genes lista
        """
        pol = {}
        mat = []
        for genis in lista:
            pol[genis] = self.get_gene(genis)[0]
        for n, m1 in enumerate(lista):
            print n
            mat.append([])
            for m2 in lista:
                ql = 0.0
                for node in self.pl:
                    if pol[m1][node] > 0.0:
                        log2m1 = numpy.log2(pol[m1][node])
                    else:
                        log2m1 = 0.0
                    if pol[m2][node] > 0.0:
                        log2m2 = numpy.log2(pol[m2][node])
                    else:
                        log2m2 = 0.0
                    if pol[m2][node]+pol[m1][node] > 0.0:
                        log2m1m2 = numpy.log2(0.5*(pol[m2][node]+pol[m1][node]))
                    else:
                        log2m1m2 = 0.0
                    ql += 0.5*pol[m1][node]*log2m1+0.5*pol[m2][node]*log2m2 - \
                          0.5*(pol[m1][node]+pol[m2][node])*log2m1m2
                mat[-1].append(numpy.sqrt(ql))
        return mat

    def adjacency_matrix(self, lista, ind=1):
        """
        Computes the Adjacency divergence matrix of the list of genes lista
        """
        cor = float(len(self.pl))/float(len(self.pl)-1)
        pol = {}
        mat = []
        for genis in lista:
            pol[genis] = self.get_gene(genis)[0]
        for n, m1 in enumerate(lista):
            print n
            ex1 = []
            for uu in self.pl:
                ex1.append(pol[m1][uu])
            ex1 = numpy.array(ex1)
            mat.append([])
            for m2 in lista:
                ex2 = []
                for uu in self.pl:
                    ex2.append(pol[m2][uu])
                ex2 = numpy.array(ex2)
                mat[-1].append(cor*numpy.dot(numpy.transpose(ex1),
                                             numpy.dot(numpy.linalg.matrix_power(self.adj, ind), ex2)))
        return mat

    def draw_net(self, color, connected=True, labels=False, ccmap='jet', weight=8.0, save='', ignore_log=False,
                 table=False):
        """
        Plots colored network
        """
        if connected:
            pg = self.gl
            pos = self.posgl
        else:
            pg = self.g
            pos = self.posg
        fig = pylab.figure()
        networkx.draw_networkx_edges(pg, pos, width=1, alpha=0.4)
        sizes = numpy.array([len(self.dic[int(node)]) for node in pg.nodes()])*weight
        values = []
        if type(color) == str or (type(color) == list and len(color) == 1):
            coloru, tol = self.get_gene(color, ignore_log=ignore_log)
            values = [coloru[node] for node in pg.nodes()]
            networkx.draw_networkx_nodes(pg, pos, node_color=values, node_size=sizes, cmap=pylab.get_cmap(ccmap))
        elif type(color) == list and len(color) == 2:
            colorr, tolr = self.get_gene(color[0], ignore_log=ignore_log)
            rmax = float(max(colorr.values()))
            if rmax == 0.0:
                rmax = 1.0
            colorb, tolb = self.get_gene(color[1], ignore_log=ignore_log)
            bmax = float(max(colorb.values()))
            if bmax == 0.0:
                bmax = 1.0
            values = [(1.0-colorb[node]/bmax, max(1.0-(colorr[node]/rmax+colorb[node]/bmax), 0.0),
                       1.0-colorr[node]/rmax) for node in pg.nodes()]
            networkx.draw_networkx_nodes(pg, pos, node_color=values, node_size=sizes)
        elif type(color) == list and len(color) == 3:
            colorr, tolr = self.get_gene(color[0], ignore_log=ignore_log)
            rmax = float(max(colorr.values()))
            if rmax == 0.0:
                rmax = 1.0
            colorg, tolg = self.get_gene(color[1], ignore_log=ignore_log)
            gmax = float(max(colorg.values()))
            if gmax == 0.0:
                gmax = 1.0
            colorb, tolb = self.get_gene(color[2], ignore_log=ignore_log)
            bmax = float(max(colorb.values()))
            if bmax == 0.0:
                bmax = 1.0
            values = [(max(1.0-(colorg[node]/gmax+colorb[node]/bmax), 0.0),
                       max(1.0-(colorr[node]/rmax+colorb[node]/bmax), 0.0),
                       max(1.0-(colorr[node]/rmax+colorg[node]/gmax), 0.0)) for node in pg.nodes()]
            networkx.draw_networkx_nodes(pg, pos, node_color=values, node_size=sizes)
        if labels:
            networkx.draw_networkx_labels(pg, pos, font_size=5, font_family='sans-serif')
        frame1 = pylab.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        if table:
            if type(color) == str:
                cell_text = [[str(min(values)*tol), str(max(values)*tol), str(numpy.median(values)*tol),
                              str(self.expr(color)), str(self.delta(color)), str(self.connectivity(color, ind=1)),
                              str(self.connectivity_pvalue(color, ind=1, c=1, n=200))]]
                columns = ['Min.', 'Max.', 'Median', 'Node fraction', 'Coeff. var.', 'Connectivity', 'p value']
                rows = [color]
                pylab.table(cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom')
            elif type(color) == list and len(color) == 1:
                cell_text = [[str(min(values)*tol), str(max(values)*tol), str(numpy.median(values)*tol),
                              str(self.expr(color[0])), str(self.delta(color[0])),
                              str(self.connectivity(color[0], ind=1)),
                              str(self.connectivity_pvalue(color[0], ind=1, c=1, n=200))]]
                columns = ['Min.', 'Max.', 'Median', 'Node fraction', 'Coeff. var.', 'Connectivity', 'p value']
                rows = [color[0]]
                pylab.table(cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom')
            elif type(color) == list and len(color) == 2:
                valuesr = colorr.values()
                valuesb = colorb.values()
                cell_text = [[str(min(valuesr)*tolr), str(max(valuesr)*tolr), str(numpy.median(valuesr)*tolr),
                              str(self.expr(color[0])), str(self.delta(color[0])),
                              str(self.connectivity(color[0], ind=1)),
                              str(self.connectivity_pvalue(color[0], ind=1, c=1, n=200))], [str(min(valuesb)*tolb),
                              str(max(valuesb)*tolb), str(numpy.median(valuesb)*tolb),
                              str(self.expr(color[1])), str(self.delta(color[1])),
                              str(self.connectivity(color[1], ind=1)),
                              str(self.connectivity_pvalue(color[1], ind=1, c=1, n=200))]]
                columns = ['Min.', 'Max.', 'Median', 'Node fraction', 'Coeff. var.', 'Connectivity', 'p value']
                rows = [color[0], color[1]]
                pylab.table(cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom', rowColours=['r', 'b'])
            elif type(color) == list and len(color) == 3:
                valuesr = colorr.values()
                valuesg = colorg.values()
                valuesb = colorb.values()
                cell_text = [[str(min(valuesr)*tolr), str(max(valuesr)*tolr), str(numpy.median(valuesr)*tolr),
                              str(self.expr(color[0])), str(self.delta(color[0])),
                              str(self.connectivity(color[0], ind=1)),
                              str(self.connectivity_pvalue(color[0], ind=1, c=1, n=200))],
                             [str(min(valuesg)*tolg), str(max(valuesg)*tolg), str(numpy.median(valuesg)*tolg),
                              str(self.expr(color[1])), str(self.delta(color[1])),
                              str(self.connectivity(color[1], ind=1)),
                              str(self.connectivity_pvalue(color[1], ind=1, c=1, n=200))],
                             [str(min(valuesb)*tolb), str(max(valuesb)*tolb), str(numpy.median(valuesb)*tolb),
                              str(self.expr(color[2])), str(self.delta(color[2])),
                              str(self.connectivity(color[2], ind=1)),
                              str(self.connectivity_pvalue(color[2], ind=1, c=1, n=200))]]
                columns = ['Min.', 'Max.', 'Median', 'Node fraction', 'Coeff. var.', 'Connectivity', 'p value']
                rows = [color[0], color[1], color[2]]
                the_table = pylab.table(cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom',
                                        rowColours=['r', 'g', 'b'], colWidths=[0.08] * 7)
                the_table.scale(1.787, 1)
                pylab.subplots_adjust(bottom=0.2)
        if save == '':
            pylab.show()
        else:
            fig.savefig(save)
        return values

    def show_statistics(self):
        """
        Shows graph statistics
        """
        f, axarr = pylab.subplots(2)
        x = map(len, self.dic.values())
        axarr[0].hist(x, max(x)-1, alpha=0.6, color='b')
        axarr[0].set_xlabel('cells per node')
        x = []
        for q in self.g.edges():
            x.append(len(set(self.dic[int(q[0])]).intersection(self.dic[int(q[1])])))
        axarr[1].hist(x, max(x)-1, alpha=0.6, color='g')
        axarr[1].set_xlabel('shared cells between connected nodes')
        pylab.show()


class RootedGraph(UnrootedGraph):
    """
    Inherits from UnrootedGraph and implements rooting methods
    """
    def get_distroot(self, root):
        """
        Dictionary of node distances to root
        """
        distroot = {}
        for i in sorted(self.pl):
            distroot[str(i)] = networkx.shortest_path_length(self.gl, str(root), i)
        return distroot

    def get_dendrite(self, permut=False):
        """
        Dendritic function
        """
        dendrite = {}
        daycolor = self.get_gene(self.rootlane, permut=permut, ignore_log=True)
        for i in self.pl:
            distroot = self.get_distroot(i)
            x = []
            y = []
            for q in distroot.keys():
                if distroot[q] != max(distroot.values()):
                    x.append(daycolor[q]-min(daycolor.values()))
                    y.append(distroot[q])
            dendrite[str(i)] = -scipy.stats.spearmanr(x, y)[0]
        return dendrite

    def find_root(self, dendritem):
        """
        Finds root and leaf
        """
        q = 1000.0
        q2 = -1000.0
        ind = 0
        ind2 = 0
        for n in dendritem.keys():
            if -2.0 < dendritem[n] < q:
                q = dendritem[n]
                ind = n
            if dendritem[n] > q2 and dendritem[n] > -2.0:
                q2 = dendritem[n]
                ind2 = n
        return ind, ind2

    def dendritic_graph(self):
        """
        Builds dendritic graph
        """
        diam = networkx.diameter(self.gl)
        g3 = networkx.Graph()
        dicdend = {}
        for n in range(diam-1):
            nodedist = []
            for k in self.pl:
                dil = networkx.shortest_path_length(self.gl, self.root, k)
                if dil == n:
                    nodedist.append(str(k))
            g2 = self.gl.subgraph(nodedist)
            dicdend[n] = sorted(networkx.connected_components(g2))
            for n2, yu in enumerate(dicdend[n]):
                g3.add_node(str(n) + '_' + str(n2))
                if n > 0:
                    for n3, yu2 in enumerate(dicdend[n-1]):
                        if networkx.is_connected(self.gl.subgraph(yu+yu2)):
                            g3.add_edge(str(n) + '_' + str(n2), str(n-1) + '_' + str(n3))
        return g3, dicdend

    def __init__(self, name, table, rootlane, shift=None, log2=False, posgl=False):
        """
        Initializes the class
        """
        UnrootedGraph.__init__(self, name, table, shift, log2, posgl)
        self.rootlane = rootlane
        self.root, self.leaf = self.find_root(self.get_dendrite())
        self.g3, self.dicdend = self.dendritic_graph()
        self.g3prun = self.g3
        self.edgesize = []
        self.edgesizeprun = []
        self.nodesize = []
        self.dicmelisa = {}
        self.nodesizeprun = []
        self.dicmelisaprun = {}
        for ee in self.g3.edges():
            yu = self.dicdend[int(ee[0].split('_')[0])][int(ee[0].split('_')[1])]
            yu2 = self.dicdend[int(ee[1].split('_')[0])][int(ee[1].split('_')[1])]
            self.edgesize.append(self.gl.subgraph(yu+yu2).number_of_edges()-self.gl.subgraph(yu).number_of_edges()
                                 - self.gl.subgraph(yu2).number_of_edges())
        for ee in self.g3.nodes():
            lisa = []
            for uu in self.dicdend[int(ee.split('_')[0])][int(ee.split('_')[1])]:
                lisa += self.dic[int(uu)]
            self.nodesize.append(len(set(lisa)))
            self.dicmelisa[ee] = set(lisa)
        try:
            self.posg3 = networkx.graphviz_layout(self.g3, 'sfdp')
            self.posg3prun = networkx.graphviz_layout(self.g3prun, 'sfdp')
        except:
            self.posg3 = networkx.spring_layout(self.g3)
            self.posg3prun = networkx.spring_layout(self.g3prun)
        self.dicdis = self.get_distroot(self.root)

    def draw_diff_net(self, color, labels=False, ccmap='jet', weight=8.0, save='', ignore_log=False):
        """
        Plots colored network
        """
        values = []
        pg = self.g3
        pos = self.posg3
        edgesize = self.edgesize
        nodesize = self.nodesize
        fig = pylab.figure()
        networkx.draw_networkx_edges(pg, pos,
                                     width=numpy.log2(numpy.array(edgesize)+1)*8.0/float(numpy.log2(1+max(edgesize))),
                                     alpha=0.6)
        if type(color) == str or (type(color) == list and len(color) == 1):
            values = []
            for _ in pg.nodes():
                values.append(0.0)
            if type(color) == str:
                color = [[color]]
            for colorm in color[0]:
                geys = self.dicgenes[colorm]
                for llp, ee in enumerate(pg.nodes()):
                    pol = 0.0
                    if self.log2 and not ignore_log:
                        for uni in self.dicmelisa[ee]:
                                pol += (numpy.power(2, float(geys[uni]))-1.0)
                        pol = numpy.log2(1.0+(pol/float(len(self.dicmelisa[ee]))))
                    else:
                        for uni in self.dicmelisa[ee]:
                            pol += geys[uni]
                        pol /= len(self.dicmelisa[ee])
                    values[llp] += pol
            networkx.draw_networkx_nodes(pg, pos, node_color=values,
                                         node_size=numpy.array(nodesize)*weight*50.0/float(max(nodesize)),
                                         cmap=pylab.get_cmap(ccmap))
        # Need to adapt next cases to average genes
        elif type(color) == list and len(color) == 2:
            geysr = self.dicgenes[color[0]]
            geysb = self.dicgenes[color[1]]
            colorr = {}
            colorb = {}
            for ee in pg.nodes():
                polr = 0.0
                polb = 0.0
                if self.log2 and not ignore_log:
                    for uni in self.dicmelisa[ee]:
                            polr += (numpy.power(2, float(geysr[uni]))-1.0)
                            polb += (numpy.power(2, float(geysb[uni]))-1.0)
                    polr = numpy.log2(1.0+(polr/float(len(self.dicmelisa[ee]))))
                    polb = numpy.log2(1.0+(polb/float(len(self.dicmelisa[ee]))))
                else:
                    for uni in self.dicmelisa[ee]:
                        polr += geysr[uni]
                        polb += geysb[uni]
                    polr /= len(self.dicmelisa[ee])
                    polb /= len(self.dicmelisa[ee])
                colorr[ee] = polr
                colorb[ee] = polb
            rmax = float(max(colorr.values()))
            bmax = float(max(colorb.values()))
            values = [(1.0-colorb[node]/bmax, max(1.0-(colorr[node]/rmax+colorb[node]/bmax), 0.0),
                       1.0-colorr[node]/rmax) for node in pg.nodes()]
            networkx.draw_networkx_nodes(pg, pos, node_color=values,
                                         node_size=numpy.array(nodesize)*weight*50.0/float(max(nodesize)))
        elif type(color) == list and len(color) == 3:
            geysr = self.dicgenes[color[0]]
            geysg = self.dicgenes[color[1]]
            geysb = self.dicgenes[color[2]]
            colorr = {}
            colorg = {}
            colorb = {}
            for ee in pg.nodes():
                polr = 0.0
                polg = 0.0
                polb = 0.0
                if self.log2 and not ignore_log:
                    for uni in self.dicmelisa[ee]:
                            polr += (numpy.power(2, float(geysr[uni]))-1.0)
                            polg += (numpy.power(2, float(geysg[uni]))-1.0)
                            polb += (numpy.power(2, float(geysb[uni]))-1.0)
                    polr = numpy.log2(1.0+(polr/float(len(self.dicmelisa[ee]))))
                    polg = numpy.log2(1.0+(polg/float(len(self.dicmelisa[ee]))))
                    polb = numpy.log2(1.0+(polb/float(len(self.dicmelisa[ee]))))
                else:
                    for uni in self.dicmelisa[ee]:
                        polr += geysr[uni]
                        polg += geysg[uni]
                        polb += geysb[uni]
                    polr /= len(self.dicmelisa[ee])
                    polg /= len(self.dicmelisa[ee])
                    polb /= len(self.dicmelisa[ee])
                colorr[ee] = polr
                colorg[ee] = polg
                colorb[ee] = polb
            rmax = float(max(colorr.values()))
            gmax = float(max(colorg.values()))
            bmax = float(max(colorb.values()))
            values = [(max(1.0-(colorg[node]/gmax+colorb[node]/bmax), 0.0),
                       max(1.0-(colorr[node]/rmax+colorb[node]/bmax), 0.0),
                       max(1.0-(colorr[node]/rmax+colorg[node]/gmax), 0.0)) for node in pg.nodes()]
            networkx.draw_networkx_nodes(pg, pos, node_color=values,
                                         node_size=numpy.array(nodesize)*weight*50.0/float(max(nodesize)))
        if labels:
            networkx.draw_networkx_labels(pg, pos, font_size=5, font_family='sans-serif')
        if save == '':
            pylab.show()
        else:
            fig.savefig(save)
        return values

    def centroid(self, genin, ignore_log=False):
        """
        Computes centroid and dispersion
        """
        dicge = self.get_gene(genin, ignore_log=ignore_log)
        pel1 = 0.0
        pel2 = 0.0
        pel3 = 0.0
        for node in self.pl:
            pel1 += self.dicdis[node]*dicge[node]
            pel2 += dicge[node]
        if pel2 > 0.0:
            cen = float(pel1)/float(pel2)
            for node in self.pl:
                pel3 += numpy.power(self.dicdis[node]-cen, 2)*dicge[node]
            return [cen, numpy.sqrt(pel3/float(pel2))]
        else:
            return [None, None]

    def get_gene(self, genin, permut=False, ignore_log=False):
        """
        Returns a dictionary that asigns to each node ID the average value of the column genin in the raw table, for
         the rows contained in that column. If permut = True performs a reshuffling of rows. If log2 = True assumes
          values in the format log2(1+x) when computing averages.
        """
        if genin == '_dist_root':
            return self.get_distroot(self.root)
        else:
            return UnrootedGraph.get_gene(self, genin, permut, ignore_log)[0]
