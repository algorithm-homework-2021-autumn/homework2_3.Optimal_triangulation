// Optimal_triangulation.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>

const double Pi = acos(-1.0);//圆周率
const double EARTH_RADIUS = 6378.137;//计算距离所需常数
const double my_eps = 1e-16;//极小值，用来判0

const int data_size = 2;//数据文件数量
const int maxn = 100;//多边形最大顶点数

class eNodeB {//基站(凸多边形顶点)类
public:
    int id;
    double lon;
    double lat;
    int order;

    eNodeB(int idx = 0, double lonx = 0, double latx = 0, int orderx = 0) {//构造函数
        id = idx;
        lon = lonx;
        lat = latx;
        order = orderx;
    }
    eNodeB(const eNodeB& b) {//拷贝构造函数
        id = b.id;
        lon = b.lon;
        lat = b.lat;
        order = b.order;
    }
    void get(int idx = 0, double lonx = 0, double latx = 0, int orderx = 0) {//赋值函数
        id = idx;
        lon = lonx;
        lat = latx;
        order = orderx;
    }
    double operator-(const eNodeB& b) {//两个基站类相减返回值为距离
        double radlon = lon * Pi / 180;
        double radlat = lat * Pi / 180;
        double radlatb = b.lat * Pi / 180;
        double radlonb = b.lon * Pi / 180;
        double v = acos(cos(radlat) * cos(radlatb) * cos(radlon - radlonb) + sin(radlat) * sin(radlatb));
        v = v * EARTH_RADIUS;
        v *= 1000;
        if (v >= 0)return v;
        else return -v;
    }
    void operator=(const eNodeB& b) {//赋值运算符重载
        id = b.id;
        lon = b.lon;
        lat = b.lat;
        order = b.order;
    }
};
void get_polygon_data(eNodeB** a, int* size); //从文件读取凸多边形数据
void getSplitP(int** s, int l, int r); //递归查找并输出每个多边形的剖分点
void get_Triangulation1(eNodeB* a, int size); //O(n^3)求凸多边形的最优三角剖分
void get_Triangulation2(eNodeB* a, int size); //O(n^2)求凸多边形的三角剖分

void get_polygon_data(eNodeB** a, int* size) {//从文件读取凸多边形数据
    std::ifstream Data;
    for (int i = 0; i < data_size; ++i) {
        std::string str = "polygon_data";
        str = str + char('1' + i);
        str = str + ".txt";
        Data.open(str);
        int j = 0;
        while(Data >> a[i][j].id) {
            Data >> a[i][j].lon;
            Data >> a[i][j].lat;
            Data >> a[i][j].order;
            ++j;
        }size[i] = j;
        Data.close();
    }
}
void getSplitP(int** s, int l, int r) {
    if (l == r)return;
    std::cout << l - 1 <<' '<< s[l][r] <<' '<< r << std::endl;
    getSplitP(s, l, s[l][r]);
    getSplitP(s, s[l][r] + 1, r);

}
void get_Triangulation1(eNodeB* a, int size) {
    std::cout << "多边形顶点数为 " << size << std::endl;
    int start_time = clock();
    double** f = new double* [size + 1]; 
    for (int i = 0; i < size + 1; ++i)f[i] = new double[size + 1];
    int** p = new int* [size + 1];
    for (int i = 0; i < size + 1; ++i)p[i] = new int[size + 1];
    //f[i][j]为凸子多边形{ vi - 1,vi,…,vj }的最优三角剖分所对应的权函数值
    //p[i][j]最优三角剖分取的分割点
    for (int i = 0; i < size; ++i)f[i][i] = 0;
    for (int r = 2; r < size; ++r) {
        for (int i = 1, j; i < size - r + 1; ++i) {
            j = i + r - 1;
            f[i][j] = -1;//初始化最小值为负数
            p[i][j] = i;
            for (int k = i; k < j; ++k) {
                double v = f[i][k]+f[k+1][j];
                v += a[i-1]-a[k];
                v += a[j]-a[k];
                v += a[i-1]-a[j];
                if (f[i][j]<0 || f[i][j]>v) {
                    f[i][j] = v;
                    p[i][j] = k;
                }
            }
        }
    }

    std::cout << "O(n^3)查找最优三角剖分的权函数值为" << std::setprecision(10) << f[1][size - 1] << std::endl;
    std::cout << "每个最优三角的三个点的序号如下:"<<std::endl;
    getSplitP(p, 1, size - 1);
    std::cout << "运行时间为 " << clock() - start_time << "ms" << std::endl;

    for (int i = 0; i < size + 1; ++i)delete[]f[i];
    delete[]f; 
    for (int i = 0; i < size + 1; ++i)delete[]p[i];
    delete[]p;
    std::cout << std::endl; 
}
void get_Triangulation2(eNodeB* a, int size) {
    std::cout << "多边形顶点数为 " << size << std::endl;
    int start_time = clock();
    double** f = new double* [size + 1];
    for (int i = 0; i < size + 1; ++i)f[i] = new double[size + 1];
    int** p = new int* [size + 1];
    for (int i = 0; i < size + 1; ++i)p[i] = new int[size + 1];
    //f[i][j]为凸子多边形{ vi - 1,vi,…,vj }的最优三角剖分所对应的权函数值
    //p[i][j]最优三角剖分取的分割点
    for (int i = 0; i < size; ++i)f[i][i] = 0;
    for (int r = 2; r < size; ++r) {
        for (int i = 1, j; i < size - r + 1; ++i) {
            j = i + r - 1;
            f[i][j] = -1;//初始化最小值为负数
            p[i][j] = i;
            int z = j - i;
            for (int w = 1; w < 4; ++w) {
                int k = i + z * w / 4;
                double v = f[i][k] + f[k + 1][j];
                v += a[i - 1] - a[k];
                v += a[j] - a[k];
                v += a[i - 1] - a[j];
                if (f[i][j]<0 || f[i][j]>v) {
                    f[i][j] = v;
                    p[i][j] = k;
                }
            }
        }
    }

    std::cout << "O(n^2)查找三角剖分的权函数值为" << std::setprecision(10) << f[1][size - 1] << std::endl;
    std::cout << "每个三角的三个点的序号如下:" << std::endl;
    getSplitP(p, 1, size - 1);
    std::cout << "运行时间为 " << clock() - start_time << "ms" << std::endl;

    for (int i = 0; i < size + 1; ++i)delete[]f[i];
    delete[]f;
    for (int i = 0; i < size + 1; ++i)delete[]p[i];
    delete[]p;
    std::cout << std::endl;
}
int main(){
    //new凸多边形数据数组和储存每组多边形点数的数组
    eNodeB**a=new eNodeB *[data_size];
    for (int i = 0; i < data_size; ++i)a[i] = new eNodeB[maxn];
    int* size = new int[data_size];
    
    //从文件读取凸多边形数据
    get_polygon_data(a,size);
    
    //求每组数据的最优三角剖分
    for (int i = 0; i < data_size; ++i) {
        std::cout << "计算第" << i + 1 << "组凸多边形的最优三角剖分如下" << std::endl;
        get_Triangulation1(a[i], size[i]);
        std::cout << "O(n^2)计算第" << i + 1 << "组凸多边形的三角剖分如下" << std::endl;
        get_Triangulation2(a[i], size[i]);
    }

    //delete分配的动态空间
    for (int i = 0; i < data_size; ++i)delete[]a[i];
    delete[]a;
    delete[]size;
    return 0;
}

