"""
Copyright 2015, Pablo G. Camara, Columbia University
"""

import networkx
import json
import numexpr
import numpy
import numpy.linalg
import numpy.random
import random
import requests
import scipy.stats
import scipy.cluster.hierarchy as sch
import scipy.interpolate
import sklearn.metrics.pairwise
import pickle
import pylab
from mpl_toolkits.mplot3d import Axes3D


"""
GLOBAL METHODS
"""


def ParseAyasdiGraph(lab, source, user, password, name, jsonout=False):
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
    if jsonout:
        with open(name + '.json', 'wb') as handle3:
            json.dump(dic2, handle3)
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
    r = session.get('https://core.ayasdi.com/v0/sources/' + source + '/networks/' + lab + '/node_groups')
    sp = json.loads(r.content)
    dicgroups = {}
    for m in sp:
        dicgroups[m['name']] = m['node_ids']
    with open(name + '.groups.pickle', 'wb') as handle:
        pickle.dump(dicgroups, handle)
    if jsonout:
        with open(name + '.json', 'wb') as handle3:
            json.dump(dicgroups, handle3)


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


def hierarchical_clustering(mat, method='median', labels=None):
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
    if labels is None:
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
    else:
        axmatrix.set_xticks(range(len(labels)))
        lab = [labels[idx1[m]] for m in range(len(labels))]
        axmatrix.set_xticklabels(lab)
        axmatrix.set_yticks(range(len(labels)))
        axmatrix.set_yticklabels(lab)
        for tick in pylab.gca().xaxis.iter_ticks():
            tick[0].label2On = False
            tick[0].label1On = True
            tick[0].label1.set_rotation('vertical')
        for tick in pylab.gca().yaxis.iter_ticks():
            tick[0].label2On = True
            tick[0].label1On = False
    axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.6])
    pylab.colorbar(im, cax=axcolor)
    pylab.show()
    return Z1


def find_clusters(z):
    clus = {}
    for y, m in enumerate(z['dcoord']):
        for n, q in enumerate(m):
            if q == 0.0:
                if z['color_list'][y] not in clus.keys():
                    clus[z['color_list'][y]] = [z['leaves'][int((z['icoord'][y][n]-5.0)/10.0)]]
                else:
                    clus[z['color_list'][y]].append(z['leaves'][int((z['icoord'][y][n]-5.0)/10.0)])
    return clus


"""
CLASSES
"""


class Preprocess(object):
    def __init__(self, fil, long, batch, cel):
        self.len = cel
        self.fil = fil
        self.long = long
        self.batch = batch
        carma = []
        self.cal = []
        self.cal2 = []
        self.tal = {}
        totalspikes = {}
        self.totaltransc = {}
        nofeature = {}
        self.spikes_ratio = {}
        for f in fil:
            totalspikes[f] = numpy.array([0.0 for _ in range(self.len)])
            self.totaltransc[f] = numpy.array([0.0 for _ in range(self.len)])
            self.spikes_ratio[f] = numpy.array([0.0 for _ in range(self.len)])
            fol = open(f, 'r')
            qty = 0
            for line in fol:
                sp = line[:-1].split('\t')
                if 'ERCC' in sp[0]:
                    q = numpy.array([float(r) for r in sp[1:]])
                    totalspikes[f] += q
                    if numpy.mean(q) > 5.0:
                        carma.append(q/numpy.mean(q))
                    qty += 1
                elif '__no_feature' in sp[0]:
                    nofeature[f] = numpy.array([float(r) for r in sp[1:]])
                elif '__ambiguous' not in sp[0] and '__alignment' not in sp[0]:
                    self.totaltransc[f] += numpy.array([float(r) for r in sp[1:]])
            self.spikes_ratio[f] = totalspikes[f]/(self.totaltransc[f] + nofeature[f])
            self.cal2 += list(self.spikes_ratio[f])
            fol.close()
            factors = []
            for k in numpy.transpose(numpy.array(carma)):
                factors.append(numpy.mean(k))
            self.cal += factors
            self.tal[f] = factors

    def show_statistics(self):
        pylab.figure()
        pylab.plot(self.cal, self.cal2, 'k.', alpha=0.6)
        pylab.yscale('log')
        pylab.xlabel('ERCC reads / average ERCC reads library')
        pylab.ylabel('ERCC reads / uniquely mapped reads')
        pylab.figure()
        q = []
        for f in self.fil:
            fol = open(f, 'r')
            for line in fol:
                sp = line[:-1].split('\t')
                for nok, tuio in enumerate(sp[1:]):
                    if float(tuio) > 0.0:
                        q.append(numpy.log2(1.0+1000000.0*float(tuio)/(self.totaltransc[f][nok])))
            fol.close()
        q = numpy.array(q)
        pylab.figure()
        pylab.hist(q, 100)
        pylab.xlabel('log_2 (1 + TPM)')
        pylab.show()

    def save(self, name, filter1low=0.0, filter1high=100000000.0, filter2low=0.0, filter2high=100000000.0,
             filter3low=0.0, filter3high=10000000000000.0):
        g = open(name, 'w')
        elim = []
        conta = 0
        cutted = {}
        for n, f in enumerate(self.fil):
            cutted[f] = 0
            genes = []
            valu = []
            valu3 = []
            for nok in range(40):
                if filter1high > self.tal[f][nok] > filter1low and filter2low < self.spikes_ratio[f][nok] < filter2high \
                        and self.totaltransc[f][nok] > 0.0:
                    valu3.append(nok)
                else:
                    elim.append(f[1])
            fol = open(f, 'r')
            for n3, line in enumerate(fol):
                sp = line[:-1].split('\t')
                if 'ERCC' not in sp[0] and '__' not in sp[0]:
                    coy = []
                    for nok, tuio in enumerate(sp[1:]):
                        if filter1high > self.tal[f][nok] > filter1low and filter2low < self.spikes_ratio[f][nok] < \
                                filter2high and self.totaltransc[f][nok] > 0.0:
                            if filter3high > float(tuio) > filter3low:
                                coy.append(numpy.log2(1.0+1000000.0*float(tuio)/(self.totaltransc[f][nok])))
                            else:
                                coy.append(0.0)
                        elif n3 == 0:
                            cutted[f] += 1
                    valu.append(numpy.array(coy))
                    genes.append(sp[0])
            fol.close()
            valu2 = numpy.transpose(numpy.array(valu))
            if n == 0:
                p = 'ID\tday\tlib\t'
                for l in genes:
                    p += l + '\t'
                g.write(p[:-1] + '\n')
            for tuj, l in enumerate(valu2):
                p = f[:-4] + '_' + str(valu3[tuj]) + '\t' + f[1] + '\t' + f[-7:-4] + '\t'
                for r in l:
                    p += str(r) + '\t'
                g.write(p[:-1] + '\n')
                conta += 1
        g.close()
        return conta, cutted


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
        with open(name + '.groups.pickle', 'r') as handler:
            self.dicgroups = pickle.load(handler)
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
        lista = map(int, self.pl)
        self.samples = []
        for i in lista:
            self.samples += self.dic[i]
        self.samples = numpy.array(list(set(self.samples)))

    def get_gene(self, genin, permut=False, ignore_log=False, con=True):
        """
        Returns a dictionary that asigns to each node ID the average value of the column genin in the raw table, for
         the rows contained in that column. If permut = True performs a reshuffling of rows. If log2 = True assumes
          values in the format log2(1+x) when computing averages.
        """
        if type(genin) != list:
            genin = [genin]
        genecolor = {}
        lista = []
        for i in self.dic.keys():
            if con:
                if str(i) in self.pl:
                    genecolor[str(i)] = 0.0
                    lista.append(i)
            else:
                genecolor[str(i)] = 0.0
                lista.append(i)
        for mju in genin:
            if mju is None:
                for i in sorted(lista):
                    genecolor[str(i)] = 0.0
            else:
                geys = self.dicgenes[mju]
                kis = range(len(geys))
                if permut:
                    random.shuffle(kis)
                for i in sorted(lista):
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

    def connectivity_pvalue(self, genin, n=500):
        """
        Returns a numpy array whose columns are p_i vectors with randomly shuffled cells
        """
        if genin is not None:
            lista = map(int, self.pl)
            jk = len(lista)
            pm = numpy.zeros((jk, n), dtype='float32')
            llm = list(numpy.arange(numpy.max(self.samples)+1)[self.samples])
            koi = {k: u for u, k in enumerate(llm)}
            geys = numpy.tile(self.dicgenes[genin], (n, 1))[:, self.samples]
            map(numpy.random.shuffle, geys)
            tot = numpy.zeros(n)
            for k, i in enumerate(lista):
                pk = geys[:, numpy.array(map(lambda x: koi[x], self.dic[i]))]
                q = pk.shape[1]
                if self.log2:
                    t1 = numexpr.evaluate('sum(2**pk - 1, 1)')/q
                    pm[k, :] = numexpr.evaluate('log1p(t1)')/0.693147
                    tot += pm[k, :]
                else:
                    pm[k, :] = numexpr.evaluate('sum(pk, 1)')/q
                    tot += pm[k, :]
            pm = pm/tot
            conn = (float(jk)/float(jk-1))*numpy.einsum('ij,ij->j', numpy.dot(self.adj, pm), pm)
            return numpy.mean(conn > self.connectivity(genin))
        else:
            return 0.0

    def connectivity(self, genis, permut=False, ind=1):
        """
        Computes order ind connectivity of column genis
        """
        dicgen = self.get_gene(genis, permut)[0]
        ex = []
        for uu in self.pl:
            ex.append(dicgen[uu])
        ex = numpy.array(ex)
        cor = float(len(self.pl))/float(len(self.pl)-1)
        return cor*numpy.dot(numpy.transpose(ex),
                             numpy.dot(numpy.linalg.matrix_power(self.adj, ind), ex))

    def delta(self, genis, group=None):
        """
        Computes coefficient of variation of column genis in the graph
        """
        per = []
        dicgen, tot = self.get_gene(genis)
        if group is not None:
            if type(group) == list:
                for k in group:
                    per += self.dicgroups[k]
                per = list(set(per))
            else:
                per = self.dicgroups[group]
            mi = [dicgen[str(node)] for node in per]
        else:
            mi = [dicgen[node] for node in self.pl]
        if numpy.mean(mi) > 0.0:
            return numpy.mean(mi)*tot, numpy.min(mi)*tot, numpy.max(mi)*tot
        else:
            return 0.0, 0.0, 0.0

    def expr(self, genis, group=None):
        """
        Computes fraction of expressed nodes for column genis
        """
        per = []
        dicgen = self.get_gene(genis)[0]
        if group is not None:
            if type(group) == list:
                for k in group:
                    per += self.dicgroups[k]
                per = list(set(per))
            else:
                per = self.dicgroups[group]
            mi = [dicgen[str(node)] for node in per]
        else:
            mi = [dicgen[node] for node in self.pl]
        return len(mi) - mi.count(0.0)

    def save(self, n=500, filtercells=0, filterexp=0.0, annotation={}):
        """
        Computes and stores all statistics in file name
        """
        pol = []
        with open(self.name + '.genes.txt', 'w') as ggg:
            cul = 'Gene\tCells\tMean\tMin\tMax\tConnectivity\tp_value\tBH p-value\t'
            for m in sorted(annotation.keys()):
                cul += m + '\t'
            ggg.write(cul[:-1] + '\n')
            lp = sorted(self.dicgenes.keys())
            for gi in lp:
                if self.expr(gi) > filtercells and self.delta(gi)[2] > filterexp:
                    pol.append(self.connectivity_pvalue(gi, n=n))
            por = benjamini_hochberg(pol)
            mj = 0
            for gi in lp:
                po = self.expr(gi)
                m1, m2, m3 = self.delta(gi)
                if po > filtercells and m3 > filterexp:
                    cul = gi + '\t' + str(po) + '\t' + str(m1) + '\t' + str(m2) + '\t' + str(m3) + '\t' + \
                          str(self.connectivity(gi)) + '\t' + str(pol[mj]) + '\t' + str(por[mj]) + '\t'
                    for m in sorted(annotation.keys()):
                        if gi in annotation[m]:
                            cul += 'Y' + '\t'
                        else:
                            cul += 'N' + '\t'
                    ggg.write(cul[:-1] + '\n')
                    mj += 1
        centr = []
        disp = []
        centr2 = []
        disp2 = []
        f = open(self.name + '.genes.txt', 'r')
        for n, line in enumerate(f):
            if n > 0:
                sp = line[:-1].split('\t')
                if float(sp[7]) <= 0.05:
                    centr.append(float(sp[1]))
                    disp.append(float(sp[5]))
                else:
                    centr2.append(float(sp[1]))
                    disp2.append(float(sp[5]))
        f.close()
        pylab.scatter(centr2, disp2, alpha=0.2, s=9, c='b')
        pylab.scatter(centr, disp, alpha=0.3, s=9, c='r')
        pylab.xlabel('cells')
        pylab.ylabel('connectivity')
        pylab.yscale('log')
        pylab.ylim(0.01, 1)
        pylab.xlim(0, max(centr+centr2))
        pylab.show()

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

    def cor_matrix(self, lista, c=1):
        """
        Computes correlation distance matrix of the list of genes lista
        """
        geys = numpy.array([self.dicgenes[mju] for mju in lista])
        return sklearn.metrics.pairwise.pairwise_distances(geys, metric='correlation', n_jobs=c)

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
                              str(self.expr(color)), str(self.connectivity(color, ind=1)),
                              str(self.connectivity_pvalue(color, n=500))]]
                columns = ['Min.', 'Max.', 'Median', 'Cells', 'Connectivity', 'p value']
                rows = [color]
                pylab.table(cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom')
            elif type(color) == list and len(color) == 1:
                cell_text = [[str(min(values)*tol), str(max(values)*tol), str(numpy.median(values)*tol),
                              str(self.expr(color[0])),
                              str(self.connectivity(color[0], ind=1)),
                              str(self.connectivity_pvalue(color[0], n=500))]]
                columns = ['Min.', 'Max.', 'Median', 'Cells', 'Connectivity', 'p value']
                rows = [color[0]]
                pylab.table(cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom')
            elif type(color) == list and len(color) == 2:
                valuesr = colorr.values()
                valuesb = colorb.values()
                cell_text = [[str(min(valuesr)*tolr), str(max(valuesr)*tolr), str(numpy.median(valuesr)*tolr),
                              str(self.expr(color[0])),
                              str(self.connectivity(color[0], ind=1)),
                              str(self.connectivity_pvalue(color[0], n=500))], [str(min(valuesb)*tolb),
                              str(max(valuesb)*tolb), str(numpy.median(valuesb)*tolb),
                              str(self.expr(color[1])),
                              str(self.connectivity(color[1], ind=1)),
                              str(self.connectivity_pvalue(color[1], n=500))]]
                columns = ['Min.', 'Max.', 'Median', 'Cells', 'Connectivity', 'p value']
                rows = [color[0], color[1]]
                pylab.table(cellText=cell_text, rowLabels=rows, colLabels=columns, loc='bottom', rowColours=['r', 'b'])
            elif type(color) == list and len(color) == 3:
                valuesr = colorr.values()
                valuesg = colorg.values()
                valuesb = colorb.values()
                cell_text = [[str(min(valuesr)*tolr), str(max(valuesr)*tolr), str(numpy.median(valuesr)*tolr),
                              str(self.expr(color[0])),
                              str(self.connectivity(color[0], ind=1)),
                              str(self.connectivity_pvalue(color[0], n=500))],
                             [str(min(valuesg)*tolg), str(max(valuesg)*tolg), str(numpy.median(valuesg)*tolg),
                              str(self.expr(color[1])),
                              str(self.connectivity(color[1], ind=1)),
                              str(self.connectivity_pvalue(color[1], n=500))],
                             [str(min(valuesb)*tolb), str(max(valuesb)*tolb), str(numpy.median(valuesb)*tolb),
                              str(self.expr(color[2])),
                              str(self.connectivity(color[2], ind=1)),
                              str(self.connectivity_pvalue(color[2], n=500))]]
                columns = ['Min.', 'Max.', 'Median', 'Cells', 'Connectivity', 'p value']
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
        f, axarr = pylab.subplots(3)
        x = map(len, self.dic.values())
        axarr[0].hist(x, max(x)-1, alpha=0.6, color='b')
        axarr[0].set_xlabel('cells per node')
        x = []
        for q in self.g.edges():
            x.append(len(set(self.dic[int(q[0])]).intersection(self.dic[int(q[1])])))
        axarr[1].hist(x, max(x)-1, alpha=0.6, color='g')
        axarr[1].set_xlabel('shared cells between connected nodes')
        pel = []
        for m in self.dic.values():
            pel += list(m)
        q = []
        for m in range(max(pel)+1):
            o = pel.count(m)
            if o > 0:
                q.append(o)
        axarr[2].hist(q, max(q), alpha=0.6, color='r')
        axarr[2].set_xlabel('number of nodes contaning the same cell')
        pylab.figure()
        r = []
        for m in self.dicgenes.keys():
            r += list(self.dicgenes[m])
        r = [k for k in r if 20 > k > 0.0]
        pylab.hist(r, 100, alpha=0.6)
        pylab.xlabel('expression')
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
        daycolor = self.get_gene(self.rootlane, permut=permut, ignore_log=True)[0]
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
        self.dicedgesize = {}
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
            self.dicedgesize[ee] = self.edgesize[-1]
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

    def select_diff_path(self):
        """
        Returns linear subgraph og the differentiation network
        """
        lista = []
        last = '0_0'
        while True:
            siz = 0
            novel = None
            for ee in self.dicedgesize.keys():
                if ((ee[0] == last and float(ee[1].split('_')[0]) > float(ee[0].split('_')[0]))
                    or (ee[1] == last and float(ee[0].split('_')[0]) > float(ee[1].split('_')[0]))) \
                        and self.dicedgesize[ee] > siz and ee not in lista:
                    novel = ee
                    siz = self.dicedgesize[ee]
            if novel is not None:
                lista.append(novel)
                if float(novel[1].split('_')[0]) > float(novel[0].split('_')[0]):
                    last = novel[1]
                else:
                    last = novel[0]
            else:
                break
        return lista

    def draw_diff_net(self, color, labels=False, ccmap='jet', weight=8.0, save='', ignore_log=False, markpath=False):
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
        if markpath:
            culer = self.select_diff_path()
            edgesize2 = [self.dicedgesize[m] for m in culer]
            networkx.draw_networkx_edges(pg, pos, edgelist=culer, edge_color='r',
                                         width=numpy.log2(numpy.array(edgesize2)+1)*8.0/float(numpy.log2(1+max(edgesize))),
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
        dicge = self.get_gene(genin, ignore_log=ignore_log)[0]
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

    def count_gene(self, genin, cond, con=True):
        """
        Return a dictionary that assigns to each node ID the fraction that takes value cond
        """
        genecolor = {}
        lista = []
        for i in self.dic.keys():
            if con:
                if str(i) in self.pl:
                    genecolor[str(i)] = 0.0
                    lista.append(i)
            else:
                genecolor[str(i)] = 0.0
                lista.append(i)
        if genin is None:
            for i in sorted(lista):
                genecolor[str(i)] = 0.0
        else:
            geys = self.dicgenes[genin]
            for i in sorted(lista):
                pol = 0.0
                for j in self.dic[i]:
                    if geys[j] == cond:
                        pol += 1.0
                genecolor[str(i)] = pol/float(len(self.dic[i]))
        tol = sum(genecolor.values())
        if tol > 0.0:
            for ll in genecolor.keys():
                genecolor[ll] = genecolor[ll]/tol
        return genecolor, tol

    def get_gene(self, genin, permut=False, ignore_log=False, con=True):
        """
        Returns a dictionary that asigns to each node ID the average value of the column genin in the raw table, for
         the rows contained in that column. If permut = True performs a reshuffling of rows. If log2 = True assumes
          values in the format log2(1+x) when computing averages.
        """
        if genin == '_dist_root':
            return self.get_distroot(self.root)
        elif genin is not None and 'time_' in genin:
            return self.count_gene(self.rootlane, float(genin.split('_')[1]))
        else:
            return UnrootedGraph.get_gene(self, genin, permut, ignore_log, con)

    def draw_expr_timeline(self, genin, ignore_log=False, path=False, noshow=False):
        """
        Plots expression of gene or genes across distance to root
        """
        distroot_inv = {}
        if not path:
            pel = self.get_distroot(self.root)
            for m in pel.keys():
                if pel[m] not in distroot_inv.keys():
                    distroot_inv[pel[m]] = [m]
                else:
                    distroot_inv[pel[m]].append(m)
        else:
            pel = self.select_diff_path()
            cali = []
            for mmn in pel:
                cali.append(mmn[0])
                cali.append(mmn[1])
            cali = list(set(cali))
            for mmn in cali:
                distroot_inv[int(mmn.split('_')[0])] = self.dicdend[int(mmn.split('_')[0])][int(mmn.split('_')[1])]
        if type(genin) != list:
            genin = [genin]
        polter = {}
        for qsd in distroot_inv.keys():
            genecolor = {}
            lista = []
            for i in self.dic.keys():
                if str(i) in distroot_inv[qsd]:
                    genecolor[str(i)] = 0.0
                    lista += list(self.dic[i])
            pol = []
            for mju in genin:
                geys = self.dicgenes[mju]
                for j in lista:
                    if self.log2 and not ignore_log:
                        pol.append(numpy.power(2, float(geys[j]))-1.0)
                    else:
                        pol.append(float(geys[j]))
            pol = map(lambda xcv: numpy.log2(1+xcv), pol)
            polter[qsd] = [numpy.mean(pol)-numpy.std(pol), numpy.mean(pol), numpy.mean(pol)+numpy.std(pol)]
        x = []
        y = []
        y1 = []
        y2 = []
        po = self.plot_rootlane_correlation(False)
        for m in sorted(polter.keys()):
            x.append((m-po[1])/po[0])
            y1.append(polter[m][0])
            y.append(polter[m][1])
            y2.append(polter[m][2])
        xnew = numpy.linspace(min(x), max(x), 300)
        ynew = scipy.interpolate.spline(x, y, xnew)
        if not noshow:
            pylab.figure()
            pylab.fill_between(xnew, 0, ynew, alpha=0.5)
            pylab.xlim(min(xnew), max(xnew))
            pylab.ylim(0.0, max(ynew)*1.2)
            pylab.xlabel(self.rootlane)
            pylab.ylabel('<log2 (1+x)>')
        return polter

    def plot_rootlane_correlation(self, doplot=True):
        """
        Plots correlation between distance to root and rootlane
        """
        pel2, tol = self.get_gene(self.rootlane, ignore_log=True)
        pel = numpy.array([pel2[m] for m in self.pl])*tol
        dr2 = self.get_distroot(self.root)
        dr = numpy.array([dr2[m] for m in self.pl])
        po = scipy.stats.linregress(pel, dr)
        if doplot:
            pylab.scatter(pel, dr, s=9.0, alpha=0.7, c='r')
            pylab.xlim(min(pel), max(pel))
            pylab.ylim(0, max(dr)+1)
            pylab.xlabel(self.rootlane)
            pylab.ylabel('distance to root node')
            xk = pylab.linspace(min(pel), max(pel), 50)
            pylab.plot(xk, po[1]+po[0]*xk, 'k--', linewidth=2.0)
            pylab.show()
        return po

    def save(self, n=500, filtercells=0, filterexp=0.0, annotation={}):
        """
        Computes and stores all statistics in file name
        """
        pol = []
        with open(self.name + '.genes.txt', 'w') as ggg:
            cul = 'Gene\tCells\tMean\tMin\tMax\tConnectivity\tp_value\tBH p-value\tCentroid\tDispersion\t'
            for m in sorted(annotation.keys()):
                cul += m + '\t'
            ggg.write(cul[:-1] + '\n')
            lp = sorted(self.dicgenes.keys())
            for gi in lp:
                if self.expr(gi) > filtercells and self.delta(gi)[2] > filterexp:
                    pol.append(self.connectivity_pvalue(gi, n=n))
            por = benjamini_hochberg(pol)
            mj = 0
            for gi in lp:
                po = self.expr(gi)
                m1, m2, m3 = self.delta(gi)
                p1, p2 = self.centroid(gi)
                if po > filtercells and m3 > filterexp:
                    cul = gi + '\t' + str(po) + '\t' + str(m1) + '\t' + str(m2) + '\t' + str(m3) + '\t' +\
                          str(self.connectivity(gi)) + '\t' + str(pol[mj]) + '\t' + str(por[mj]) +\
                          '\t' + str(p1) + '\t' + str(p2) + '\t'
                    for m in sorted(annotation.keys()):
                        if gi in annotation[m]:
                            cul += 'Y' + '\t'
                        else:
                            cul += 'N' + '\t'
                    ggg.write(cul[:-1] + '\n')
                    mj += 1
        centr = []
        disp = []
        centr2 = []
        disp2 = []
        f = open(self.name + '.genes.txt', 'r')
        for n, line in enumerate(f):
            if n > 0:
                sp = line[:-1].split('\t')
                if float(sp[7]) <= 0.05:
                    centr.append(float(sp[1]))
                    disp.append(float(sp[5]))
                else:
                    centr2.append(float(sp[1]))
                    disp2.append(float(sp[5]))
        f.close()
        pylab.scatter(centr2, disp2, alpha=0.2, s=9, c='b')
        pylab.scatter(centr, disp, alpha=0.3, s=9, c='r')
        pylab.xlabel('cells')
        pylab.ylabel('connectivity')
        pylab.yscale('log')
        pylab.ylim(0.01, 1)
        pylab.xlim(0, max(centr+centr2))
        fig = pylab.figure()
        ax = fig.add_subplot(111, projection='3d')
        f = open(self.name + '.genes.txt', 'r')
        for n, line in enumerate(f):
            if n > 0:
                sp = line[:-1].split('\t')
                ax.scatter(float(sp[8]), float(sp[9]), float(sp[1]), c='k', alpha=0.2, s=10)
        ax.set_xlabel('Centroid')
        ax.set_ylabel('Dispersion')
        ax.set_zlabel('Cells')
        pylab.show()

    def candidate_subpopulations(self, lista, thres=0.05):
        mat = self.adjacency_matrix(lista)
        mat2 = self.JSD_matrix(lista)
        x = []
        y = []
        dic = {}
        n0 = 0
        for n1, m1 in enumerate(mat2):
            for n2, m2 in enumerate(m1[n1+1:]):
                pel = mat[n1][n1+n2+1]
                if m2 > 0.6 and pel > 0.05:
                    y.append(pel)
                    x.append(m2)
                    dic[n0] = (n1, n1+n2+1)
                    n0 += 1
        xr = []
        yr = []
        xb = []
        yb = []
        pairs = []
        pvalues = []
        lista2 = []
        for n, m in enumerate(x):
            count = 0.0
            for nn, t in enumerate(x):
                if t > m and y[nn] > y[n]:
                    count += 1.0
            count = float(count)/float(len(x))
            p1 = lista[dic[n][0]]
            p2 = lista[dic[n][1]]
            pairs.append((p1, p2))
            pvalues.append(count)
        pvaluescor = benjamini_hochberg(pvalues)
        for n, m in enumerate(pvaluescor):
            if m <= thres:
                xr.append(x[n])
                yr.append(y[n])
                lista2.append(pairs[n][0])
                lista2.append(pairs[n][1])
            else:
                xb.append(x[n])
                yb.append(y[n])
        lista2 = list(set(lista2))
        pylab.scatter(xb, yb, alpha=0.1, s=6.0, c='b')
        pylab.scatter(xr, yr, alpha=0.4, s=6.0, c='r')
        pylab.xlim(0.6, 1)
        pylab.ylim(0.05, 1)
        pylab.xlabel('Jensen-Shannon distance')
        pylab.ylabel('Adjacency')
        mat3 = []
        dicmat3 = []
        for n, m1 in enumerate(lista2):
            mat3.append([])
            dicmat3.append(m1)
            for m2 in lista2:
                mat3[-1].append(mat2[lista.index(m1)][lista.index(m2)])
        return find_clusters(hierarchical_clustering(mat3, labels=dicmat3, method='centroid'))
