#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <string.h>
using namespace std;
// 节点个数
//             "model1":1039722,
//             "model2":1083568,
//             "model3":997663,
//             "model4":1001733,
//             "real1":1694616,
//             "real2":1957027,
//             "real3":426485,
//             "real4":855802

// 边条数
//             "model1":4158873,
//             "model2":3250524,
//             "model3":4988570,
//             "model4":3005190,
//             "real1":11094209,
//             "real2":2760388,
//             "real3":8543321,
//             "real4":4291352

string fnames[8] = {"networks/model1.csv","networks/model2.csv",
                    "networks/model3.csv","networks/model4.csv",
                    "networks/real1.csv","networks/real2.csv",
                    "networks/real3.csv","networks/real4.csv"};
int nodenums[8] = {1039722,1083568,997663,1001733,1694616,1957027,426485,855802};
int edgenums[8] = {4158873,3250524,4988570,3005190,11094209,2760388,8543321,4291352};
int **nodeseqs = (int**) malloc(sizeof(int*) * 8); // 8个网络的节点重要性序列

//初始化clusters
int initClusters(int** clusters,int** freelist, int nodenum)
{
    int findex = 0;//freelist索引号
    int d = 1;
    int index = 0;
    while(nodenum/d>8)
    {
        //把每一个cluster的大小再翻倍 避免cluster数量不足
        int csize = nodenum/d * 2;
        //freelist每一层结构 ---> [0存放块大小 1存放尾指针 2开始存放空的块号]
        freelist[findex] = (int*) malloc(sizeof(int)*(2+d-index+1));
        freelist[findex][0] = csize;
        freelist[findex][1] = d - index + 2;
        for(int i=index;i<=d;i++)
        {
            //给cluster申请内存
            clusters[i] = (int*) malloc(sizeof(int)*(2+csize));
            //0实际大小 1逻辑大小
            clusters[i][0] = csize;
            clusters[i][1] = 0;
            //将该cluster的索引号添加到freelist上
            freelist[findex][i-index+2] = i;
        }
        findex++;
        index = d+1;
        d = d*2;
    }

    //剩下的都是规模比较小的集团 数量也比较多 直接指定大小为8
    freelist[findex] = (int*) malloc(sizeof(int)*(2+nodenum-index));
    freelist[findex][0] = 8;
    freelist[findex][1] = nodenum - index + 1;
    for(int i=index;i<nodenum;i++)
    {
        clusters[i] = (int*) malloc(sizeof(int)*(8+2));
        clusters[i][0] = 8;
        clusters[i][1] = 0;
        freelist[findex][i-index+2] = i;
    }

    //freelist逆序 把小集团放在前面
    for(int i=0;i<(findex+1)/2;i++)
    {
       int* tmp = freelist[i];
       freelist[i] = freelist[findex-i];
       freelist[findex-i] = tmp;
    }

    //freelist结束标志
    findex++;
    freelist[findex] = (int*) malloc(sizeof(int)*2);
    freelist[findex][0] = 0;
    freelist[findex][1] = 0;
    return 0;
}

//释放一个cluster 并非真的使用free()释放 只是把该cluster的ID号放回freelist列表
int freeCluster(int** clusters, int** freelist, int clusterid)
{
    //逻辑大小重置为0
    clusters[clusterid][1] = 0;
    //寻找插入位置
    int csize = clusters[clusterid][0];
    int findex = 0;
    while(csize != freelist[findex][0]) findex++;
    //修改尾指针 并插入该cluster的ID号
    int tail = freelist[findex][1];
    tail++;
    freelist[findex][1] = tail;
    freelist[findex][tail] = clusterid; 
}

//根据申请的大小 从freelist中找出一个符合条件的cluster 返回其ID号
int getFreeClusterID(int** clusters, int** freelist, int requiresize)
{
    //定位到大小合适的索引
    int findex = 0;
    while(requiresize > freelist[findex][0]) findex++;
    //从起始点开始 最多找3层(可能存在当前层级所有cluster均分配完毕 此时可以从下一层分配一个更大的)
    int s = findex;
    int sl = 3;
    while(sl>0 && freelist[s][0]>0)
    {
        int tail = freelist[s][1];
        if(tail > 1) //存在空的cluster 将其卸下
        {
            int r = freelist[s][tail];
            tail--;
            freelist[s][1] = tail;
            return r;
        }
        s++;
        sl--;
    }
    //无法分配cluster 返回错误
    cout<<"get free cluster error!\n";
    exit(-1);
}

//动态拓展cluster大小(cluster大小不够的情况下 重新分配一个更大的cluster)
//函数只返回新的cluster的ID号 并把旧的clusterid拷贝过来 最后释放旧的cluster
//node与clusterID的对应关系需要另外再更新
int extendCluster(int** clusters, int** freelist, int tid, int requiresize)
{
    int* tcluster = clusters[tid];
    //找一个大小合适的cluster
    int nid = getFreeClusterID(clusters, freelist, requiresize);
    //将原cluster的内容拷贝过来 第一个元素不用拷贝(表示块的实际大小)
    memcpy(clusters[nid]+1,tcluster+1,sizeof(int)*(tcluster[1]+1));
    //释放原cluster
    freeCluster(clusters,freelist,tid);
    return nid;
}

int insertNode(int nodeid, int maxclustersize, int** clusters, int* node_clusterid, int** freelist, int** adjmatrix)
{
    int* neibs = adjmatrix[nodeid];//得到所有邻居
    if(neibs[0]>0)
    {
        int* tcluster;
        int tid = -1;
        bool newcluster=true;//表示是否需要新建集团

        for(int j=1;j<=neibs[0];j++)//遍历邻居 找出最大集团
        {
            int pid = node_clusterid[neibs[j]];
            if(pid > -1)//若邻居在某个集团中
            {
                if(tid == -1)
                {
                    tid = pid;
                    tcluster = clusters[tid];
                }
                else//说明不止一个邻居有所属集团
                {
                    int* pcluster = clusters[pid];
                    if(tcluster[1]<pcluster[1])//尽量找大的集团 把小集团往大集团里合并
                    {
                        tid = pid;
                        tcluster = pcluster;
                    }
                }
            }
        }


        if(tid > -1)//不需要新建集团 可以将当前点加入该集团
        {
            newcluster = false;

            //把当前点添加到集团中
            if(tcluster[1]+1 > tcluster[0])//判断当前空间是否足够
            {
                //如果不够需要申请一个更大的集团
                tid = extendCluster(clusters,freelist,tid,tcluster[1]+1);
                tcluster = clusters[tid];
                //修改节点与集团ID对应关系
                for(int k=2;k<=tcluster[1]+1;k++) node_clusterid[tcluster[k]] = tid;
            }
            tcluster[1] = tcluster[1] + 1;
            tcluster[tcluster[1]+1] = nodeid;
            node_clusterid[nodeid] = tid;


            //合并邻居中的集团
            for(int j=1;j<=neibs[0];j++)
            {
                int pid = node_clusterid[neibs[j]];
                if(pid>-1 && pid != tid)//避免重复合并
                {
                    int* pcluster = clusters[pid];

                    //判断是否需要扩展大小
                    if(tcluster[1]+pcluster[1] > tcluster[0])//判断当前空间是否足够
                    {
                        //如果不够需要申请一个更大的集团
                        tid = extendCluster(clusters,freelist,tid,tcluster[1]+pcluster[1]);
                        tcluster = clusters[tid];
                        //修改节点与集团ID对应关系
                        for(int k=2;k<=tcluster[1]+1;k++) node_clusterid[tcluster[k]] = tid;
                    }

                    //合并&&修改邻居集团的集团号
                    for(int k=2;k<=pcluster[1]+1;k++)
                    {
                        tcluster[tcluster[1]+k] = pcluster[k];
                        node_clusterid[pcluster[k]] = tid;
                    }
                    tcluster[1] = tcluster[1] + pcluster[1];
                    
                    //释放原来的集团(实际就是放回freelist)
                    freeCluster(clusters,freelist,pid);
                }
            }

            clusters[tid] = tcluster;
            //更新最大集团规模
            if(maxclustersize < tcluster[1]) maxclustersize = tcluster[1];
        }

        if(newcluster)//需要新建集团
        {
            int nid = getFreeClusterID(clusters,freelist,1);
            clusters[nid][1] = 1;
            clusters[nid][2] = nodeid;
            node_clusterid[nodeid] = nid;

            if(maxclustersize < 1) maxclustersize = 1;
        }
    }
    else
    {
        cout<<"isolate node.\n";
    }
    return maxclustersize;
}

//计算鲁棒值
float calculateRobustness(int nodenum, int* nodeseq, int** adjmatrix)
{
    //存储所有集团
    //clusters[集团ID号] --> 集团[0:实际大小 1:逻辑大小 2->n:节点ID号]
    int** clusters = (int**) malloc(sizeof(int*)*(nodenum));

    //效率和内存消耗折中 把所有内存申请好 这样在后续的计算过程中不用再使用malloc,realloc折腾内存了
    //实现上模仿了Linux堆内存管理
    //存储所有空集团列表 存在freelist中
    //freelist:
    //第i层 --> [0:本层cluster块大小 1:尾指针 2->n:空集团ID号] 尾指针指向最后一个空集团ID号
    int** freelist = (int**) malloc(sizeof(int*)*30);
    //初始化
    initClusters(clusters, freelist, nodenum);

    //存储所有节点和集团ID的对应关系 初始值全为-1
    //节点ID --> 集团ID
    int* node_clusterid = (int*) malloc(sizeof(int)*nodenum);
    memset(node_clusterid, -1, sizeof(int)*nodenum);

    long sumall = 0;//总sum
    int maxclustersize = 0;//当前最大集团规模

    //反向遍历输入序列
    for(int i=nodenum-1;i>=0;i--)
    {
        maxclustersize = insertNode(nodeseq[i],maxclustersize,clusters,node_clusterid,freelist,adjmatrix);
        sumall += maxclustersize;
    }

    //释放内存
    for(int i=0;i<nodenum;i++)
    {
        free(adjmatrix[i]);
    }
    free(adjmatrix);
    for(int i=0;i<nodenum;i++)
    {
        if(clusters[i] != NULL)
        {
            free(clusters[i]);
            clusters[i] = NULL;
        }
    }
    free(clusters);
    free(node_clusterid);

    //输出鲁棒值
    cout<<sumall<<"\n";
    double auc,sump;
    sump = 1.0*(sumall-nodenum);
    auc = sump/nodenum/nodenum;
    cout << setprecision(10) << auc << "\n";
    return auc;
}


int readAllSequence(string filename)
{
    fstream infile(filename.c_str());
    string line;
    int linenum = 0;
    while(getline(infile,line))
    {
        nodeseqs[linenum] = (int*) malloc(sizeof(int) * nodenums[linenum]);
        stringstream linestream(line);
        string value;
        int colnum = 0;
        while(getline(linestream,value,','))
        {
            colnum++;
            if(colnum == 1) continue;
            int v = atoi(value.c_str());
            nodeseqs[linenum][colnum-2] = v;
        }
        linenum++;
    }
    return 0;
}

int** readGraph(string filename, int nodenum, int edgenum)
{
    int *source = (int*) malloc(sizeof(int) * edgenum);
    int *target = (int*) malloc(sizeof(int) * edgenum);
    int **adjmatrix = (int**) malloc(sizeof(int*) * nodenum); //adjmatrix 指针数组 malloc(20) 第一个元素存个数
    int EX_SIZE = 500;

    fstream infile(filename.c_str());
    int a,b;
    char comma;
    int count=0;

    while(infile >> a >> comma >> b)
    {
        source[count] = a;
        target[count] = b;
        count++;
    }

    for(int i=0;i<nodenum;i++)
    {
        adjmatrix[i] = (int*) malloc(sizeof(int) * EX_SIZE);
        adjmatrix[i][0] = 0;
    }

    for(int i=0;i<edgenum;i++)
    {
        int* a = adjmatrix[source[i]];
        int* b = adjmatrix[target[i]];

        if(a[0]%EX_SIZE == EX_SIZE-1)
        {
            a = (int*) realloc(a,sizeof(int) * (a[0]+1+EX_SIZE));
        }
        if(b[0]%EX_SIZE == EX_SIZE-1)
        {
            b = (int*) realloc(b,sizeof(int) * (b[0]+1+EX_SIZE));
        }

        //添加邻接表关系 修改第一个元素个数
        int anum = a[0] + 1;
        int bnum = b[0] + 1;
        a[0] = anum;
        b[0] = bnum;
        a[anum] = target[i];
        b[bnum] = source[i];

        adjmatrix[source[i]] = a;
        adjmatrix[target[i]] = b;
    }

    //释放内存
    free(source);
    free(target);
    return adjmatrix;
}

void printAdjmatrix(int** adjmatrix)
{
    for(int i=0;i<426485;i++)
    {
        int* p= adjmatrix[i];
        for(int j=1;j<=p[0];j++)
        {
            cout<<p[j]<<" ";
        }
        cout<<"\n";
    }
}

int main(int argc, char **argv)
{
    char *fcsv = argv[1];
    readAllSequence(fcsv);
    cout<<"Read all sequence done.\n";

    double finalresult;
    for(int i=0;i<8;i++)
    {
        cout<<"read Graph: "<<fnames[i]<<"\n";
        int** dj = readGraph(fnames[i],nodenums[i],edgenums[i]);
        cout<<"Computing...\n";
        finalresult += calculateRobustness(nodenums[i],nodeseqs[i],dj);
    }
    cout<<"final result: "<<finalresult<<"\n";
    //释放内存
    for(int i=0;i<8;i++)
    {
        free(nodeseqs[i]);
    }
    free(nodeseqs);
    return 0;
}
