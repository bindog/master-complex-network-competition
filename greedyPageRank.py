from igraph import Graph
import numpy as np
import sys
import argparse

databasepath = "networks/"
networkname = ["model1","model2","model3","model4","real1","real2","real3","real4"]

def convertTrueID(g,pr_index):
    t_index = []
    for i in pr_index:
        t_index.append(g.vs[i]["name"])
    return t_index

def greedyPageRank(g,batch_size=500):
    finallist = []
    count = 0
    while len(g.vs) > 0:
        print("Iteration:" + str(count) + "...")
        count += 1
        prlist = g.pagerank(damping=0.9)
        pr_index = np.argsort(prlist)
        pr_index = pr_index[::-1]
        if len(pr_index) <= batch_size:
            print("over")
            finallist.extend(convertTrueID(g,pr_index))
            g.delete_vertices(pr_index)
        else:
            finallist.extend(convertTrueID(g,pr_index[:batch_size]))
            g.delete_vertices(pr_index[:batch_size])
    return finallist


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Greedy PageRank", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-n", "--networkid", required=True, type=int,
    help="0 --> model1\t4 --> real1\n1 --> model2\t5 --> real2\n2 --> model3\t6 --> real3\n3 --> model4\t7 --> real4\n")
    parser.add_argument("-o", "--outfile", required=True, type=str)
    parser.add_argument("-b", "--batchsize", default=500, type=int)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    with open(args.outfile,"w") as outfile:
        fm = networkname[args.networkid]
        rawfile = open(databasepath + fm + '.csv','r')
        els = rawfile.read().replace(',', ' ')
        nfile = open(fm + '.csv.space', 'w')
        nfile.write(els)
        rawfile.close()
        nfile.close()
        print("Formatting Done.")

        elfile = open(fm + '.csv.space', 'r')
        print("Reading Graph : " + fm)
        rgraph = Graph.Read_Edgelist(elfile, directed=False)
        elfile.close()
        for i in range(len(rgraph.vs)):
            rgraph.vs[i]["name"] = i

        print("Greedy PageRank...")
        s_index = greedyPageRank(rgraph,args.batchsize)
        fstr = fm + "," + ",".join(map(str,s_index)) + "\n"
        outfile.write(fstr)
