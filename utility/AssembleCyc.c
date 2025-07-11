#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h> // 为 isspace() 函数添加的头文件

// 调试函数，确保在程序崩溃前所有信息都被打印出来
void debug_flush() {
    fflush(stdout);
}

// *** 新增的辅助函数：删除字符串末尾的空白字符 ***
void trim_trailing_whitespace(char *str) {
    if (str == NULL) return;
    int len = strlen(str);
    // 从后往前遍历，只要是空白字符，就将其替换为字符串结束符 '\0'
    while (len > 0 && isspace((unsigned char)str[len - 1])) {
        len--;
    }
    str[len] = '\0';
}


typedef struct{
    int start;
    int N;  
    int CA;
    int C;
    int O;
    int sstart;
    int end;
    char nam[10]; // 足够存储9个字符的残基名 + '\0'
}RES;

#define MAX_RES_IN_FRAG 50 
typedef struct{
    int an;
    float xyz[100][3]; // 最大100个原子
    char nam[100][4];  // 原子名，例如 "N  ", "CA ", "O  " + '\0'
    int resn;
    RES res[MAX_RES_IN_FRAG]; 
}FRAGMENT;


#define MAX_RES_IN_PEP 150 // 最大150个残基
typedef struct{
    float N[3]; 
    float CA[3];
    float C[3];
    float O[3];
    float H[3];
    char resnam[10]; // 例如 "ALA A   1" + '\0'
    int satn; // 侧链原子数量
    float side[10][3]; // 最大10个侧链原子
    char nam[10][4]; // 侧链原子名 + '\0'
}CYCRES;

#define MAXHBN  15 // 最大15个氢键
typedef struct{
    char linkname[2][100];
    int typ[2];
    float linkscore[2];
    int resn;
    int len[4];
    int hbn;
    int index[MAXHBN][2];
    float hbdis[MAXHBN];
    CYCRES res[MAX_RES_IN_PEP]; // 存储肽链的残基信息
}CYCPEP;

/*****************几何计算*************************/
void getdx(float x1[3],float x2[3],float dx[3])
{
    dx[0]=x1[0]-x2[0];
    dx[1]=x1[1]-x2[1];
    dx[2]=x1[2]-x2[2];
}
float get_len(float x[3])
{
    float r;

    r=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

    return sqrt(r);
}

void direction( float a[3], float b[3], float c[3])
{
    int i;
    float r;
    
    getdx(a,b,c);
    r=get_len(c);
    if (r < 1e-6) return; // 避免除以零或非常小的数
    for(i=0;i<3;i++) 
        c[i]=c[i]/r;
}

#define PI 3.14159265358979323846
#define TODEG(A)     ((A)*57.295779513)

float iprod(float a[3],float b[3])
{
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}
float cos_angle(const float a[3],const float b[3])
{
  float   cos_val; // 避免与标准库的cos函数名冲突
  int    m;
  double aa,bb,ip,ipa,ipb; 
  
  ip=ipa=ipb=0.0;
  for(m=0; m<3; m++) {
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  if (ipa < 1e-9 || ipb < 1e-9) return 1.0; // 避免除以零，返回1.0表示向量共线
  cos_val=ip/sqrt(ipa*ipb);
  if (cos_val > 1.0) 
    return  1.0; 
  if (cos_val <-1.0) 
    return -1.0;
  
  return cos_val;
}
float cal_angle(float x1[3], float x2[3], float x3[3])
{
    float x12[3], x32[3];
    float cosa;
    float a;

    getdx(x1,x2,x12);
    getdx(x3,x2,x32);
    cosa=cos_angle(x12,x32);
    a=acos(cosa);

    return TODEG(a);
}
void oprod(const float a[3],const float b[3],float c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

float cal_dih(float x1[3], float x2[3], float x3[3], float x4[3])
{
    float x12[3],x32[3],x34[3], m[3], n[3];
    float ipr,phi,cos_phi,sign;

    getdx(x1,x2,x12);  
    getdx(x3,x2,x32);   
    getdx(x3,x4,x34);   

    oprod(x12,x32,m); 
    oprod(x32,x34,n);
    cos_phi=cos_angle(m,n);
    phi=acos(cos_phi);
    ipr=iprod(x12,n);
    sign=(ipr<0.0)?-1.0:1.0;
    phi=sign*phi; 

    return TODEG(phi);
}

#define HBONDLEN 0.98
void calHNxyz(float H[3], float N[3], float CA[3], float C[3])
{
    float v1[3],v2[3],v3[3];
    int i;
    float len;

    direction(CA, N, v1);
    direction(C, N, v2);
    for(i=0;i<3;i++){
        v3[i]=v1[i]+v2[i];
    }
    len=get_len(v3);
    if (len < 1e-6) { // 如果v3接近零向量，则给定一个默认方向以避免除以零
        H[0] = N[0] - HBONDLEN; H[1] = N[1]; H[2] = N[2];
        return;
    }
    for(i=0;i<3;i++){
        v3[i]=v3[i]/len;
    }
    for(i=0;i<3;i++) H[i]=N[i]-v3[i]*HBONDLEN;
}
void addHN(CYCPEP *p)
{
    int ir;
    if (p->resn <= 0) {
        return; // 无残基可处理
    }
    if (p->resn > MAX_RES_IN_PEP) {
        printf("致命错误 (addHN): 肽链有 %d 个残基，但允许的最大值为 %d。\n", p->resn, MAX_RES_IN_PEP);
        debug_flush();
        exit(1);
    }

    // 假设第一个残基的N原子连接着上一个残基的C原子，如果需要处理循环肽链的特殊情况
    // 这里使用了 p->res[p->resn-1].C 作为第一个残基N的H计算的C原子
    calHNxyz(p->res[0].H, p->res[0].N, p->res[0].CA, p->res[p->resn-1].C); 
    for(ir=1;ir<p->resn;ir++){
        calHNxyz(p->res[ir].H, p->res[ir].N, p->res[ir].CA, p->res[ir-1].C);
    }
}

float distance2(float x1[3], float x2[3])
{
        return (x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2]);
}


/***********************氢键计算*******************************************/
#define HBDISMAX  (2.5*2.5) // 氢键距离上限平方
#define HBDISMIN  (1.6*1.6) // 氢键距离下限平方

int dispair(CYCPEP *p, int pair[][2], float *hbdis)
{
    int ir,jr;
    float dis;
    int n=0; // 找到的可能氢键数
    
    for(ir=0;ir<p->resn;ir++){
        for(jr=0;jr<p->resn;jr++){
            if(jr==ir) continue; // 不与自身形成氢键
            dis=distance2(p->res[ir].H, p->res[jr].O);
            if(dis<HBDISMAX&&dis>HBDISMIN){
                if (n >= MAXHBN) {
                    printf("警告 (dispair): 找到的可能氢键数超过最大值 MAXHBN (%d)，多余的将被忽略。\n", MAXHBN);
                    debug_flush();
                    return n; // 返回当前已找到的有效数量
                }
                pair[n][0]=ir;
                pair[n][1]=jr;
                hbdis[n]=sqrt(dis);
                n++;
            }
        }
    }
    return n;
}
int checkangle(CYCPEP *p, int n, int pair[][2], int *index)
{
    float ang1,ang2;
    float dih;
    int i;
    int ir,jr;
    int c=0; // 有效氢键计数
    
    for(i=0;i<n;i++){
        ir=pair[i][0];
        jr=pair[i][1];
        // 角度检查：N-H...O 角度 (氢键供体角度)
        ang1=cal_angle(p->res[ir].N, p->res[ir].H, p->res[jr].O);
        if(ang1<90) continue; // 通常要求供体角度大于90度
        // 角度检查：H...O=C 角度 (氢键受体角度)
        ang2=cal_angle(p->res[ir].H, p->res[jr].O, p->res[jr].C);
        if(ang2<90) continue; // 通常要求受体角度大于90度
        // 二面角检查：H-O-C-CA 二面角
        dih=cal_dih(p->res[ir].H, p->res[jr].O, p->res[jr].C, p->res[jr].CA);
        if(dih>50&&dih<130) continue; // 通常要求二面角在特定范围内

        if (c >= MAXHBN) {
            printf("警告 (checkangle): 找到的有效氢键数超过最大值 MAXHBN (%d)，多余的将被忽略。\n", MAXHBN);
            debug_flush();
            return c; // 返回当前已找到的有效数量
        }
        index[c]=i; // 记录有效氢键的原始索引
        c++;
    }
    return c;   
}
/*************************************************************************/

void getxyz(char line[], float xyz[3])
{
    int i;
    for(i=0;i<3;i++)
        xyz[i]=atof(line+30+8*i);
}

// *** 修改：增加行缓冲区大小以容纳长文件路径 ***
#define LINELEN 512 

void readfrag(char *pn, FRAGMENT *p)
{
    FILE *pf;
    char line[LINELEN];
    int n=0; // 原子计数 (在FRAGMENT内的索引)
    int rn=0; // 残基计数 (在FRAGMENT内的索引)
    
    if((pf=fopen(pn,"r"))==NULL){
        printf("\n致命错误: 无法打开片段结构文件 '%s'。请检查路径和文件权限。\n\n",pn);
        debug_flush();
        exit(1);
    }
    
    while(fgets(line, LINELEN,pf)){
        if(!strncmp(line,"ATOM",4)){
            // 检测新残基的逻辑
            // 如果是第一个残基，或者当前残基名与上一个不同，则视为新残基
            if(rn == 0 || strncmp(line+17,p->res[rn-1].nam,9) != 0){ 
                 if (rn >= MAX_RES_IN_FRAG) {
                    printf("致命错误 (readfrag): 文件 '%s' 中的残基数过多 (%d)，最大允许 %d。\n", pn, rn, MAX_RES_IN_FRAG);
                    debug_flush();
                    fclose(pf);
                    exit(1);
                }
                strncpy(p->res[rn].nam,line+17,9);
                p->res[rn].nam[9]='\0'; // 确保字符串以空字符结尾
                p->res[rn].start=n;
                p->res[rn].sstart=n+4; // 假设侧链原子通常在主链N,CA,C,O之后
                if(rn!=0) p->res[rn-1].end=n-1; // 更新前一个残基的结束原子索引
                rn++;
            }

            // 原子数量检查
            if (n >= 100) { // FRAGMENT的xyz和nam数组大小是100
                 printf("致命错误 (readfrag): 文件 '%s' 中的原子数过多 (%d)，最大允许 100。\n", pn, n);
                 debug_flush();
                 fclose(pf);
                 exit(1);
            }
            getxyz(line,p->xyz[n]);
            strncpy(p->nam[n],line+13,3);
            p->nam[n][3]='\0'; // 确保字符串以空字符结尾
            
            int current_res_index = rn - 1;
            if (current_res_index < 0) { // 安全检查，确保rn至少为1
                printf("警告 (readfrag): 在解析ATOMS时current_res_index为负值。文件: %s, 行: %s", pn, line);
                debug_flush();
                continue; 
            }
            if(!strncmp(line+13,"N  ",3)){
                p->res[current_res_index].N=n;
            }
            else if(!strncmp(line+13,"CA ",3)){
                p->res[current_res_index].CA=n;
            }
            else if(!strncmp(line+13,"C  ",3)){
                p->res[current_res_index].C=n;
            }
            else if(!strncmp(line+13,"O  ",3)){
                p->res[current_res_index].O=n;
            }
            else if(!strncmp(line+13,"CB ",3)){
                p->res[current_res_index].sstart=n; // 更新侧链起始索引，CB通常是第一个侧链原子
            }
            n++;
        }
    }
    fclose(pf);
    p->an=n; // 总原子数
    p->resn=rn; // 总残基数
    if (rn > 0) p->res[rn-1].end=n-1; // 更新最后一个残基的结束原子索引
}

// 辅助函数：从 PDB REMARK 行中提取文件名
void getfname(char *line, char *linkname)
{
    // 假设文件名在 REMARK link 后面，并以空格分隔
    char *start = strstr(line, "link");
    if (start) {
        start += strlen("link");
        // 跳过空格
        while (*start == ' ' && *start != '\0') {
            start++;
        }
        
        // *** 修改：在复制前，先清理字符串末尾的空白字符 ***
        trim_trailing_whitespace(start);

        strncpy(linkname, start, 99);
        linkname[99] = '\0';

    } else {
        linkname[0] = '\0'; // 未找到 "link" 关键字
    }
}

#define BESTCLUSTERNUM 3
void readlink(char *pn, FRAGMENT p[2][BESTCLUSTERNUM], int typ[2][BESTCLUSTERNUM], char fragname[2][100], char linkname[2][BESTCLUSTERNUM][100], float fragscore[2], float linkscore[2][BESTCLUSTERNUM], int tertyp[2][BESTCLUSTERNUM][2],int clusterN[2])
{
    FILE *pf;
    char line[LINELEN];
    int n=0; // 原子计数 (在当前片段内的索引)
    int rn=0; // 残基计数 (在当前片段内的索引)
    int id = 0; // 片段ID (0 或 1)
    int lid = 0; // 连接器ID (0 到 BESTCLUSTERNUM-1)
    int hasgetfragname1=0; // 标记是否已获取fragname[0]
    int hasgetfragname2=0; // 标记是否已获取fragname[1]
    
    clusterN[0] = 0; clusterN[1] = 0; // 初始化簇数量

    if((pf=fopen(pn,"r"))==NULL){
        printf("\n致命错误: 无法打开蛋白质链接文件 '%s'。请检查路径和文件权限。\n\n",pn);
        debug_flush();
        exit(1);
    }

    while(fgets(line, LINELEN,pf)){
        if(!strncmp(line,"TITLE",5)){
            // 解析 TITLE 行：TITLE ID LID Type
            // 例如 "TITLE    1 1 0" -> id=0, lid=0, typ=0
            if (sscanf(line + 8, "%d %d", &id, &lid) >= 2) { 
                id = id - 1; // 转换为0-based索引
                lid = lid - 1; // 转换为0-based索引
                
                if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                    // 安全地解析类型
                    char* type_str_start = line + 13;
                    while(*type_str_start == ' ' && *type_str_start != '\0') type_str_start++;
                    typ[id][lid] = atoi(type_str_start);
                    
                    if (clusterN[id] < lid + 1) clusterN[id] = lid + 1; // 更新簇数量
                } else {
                    printf("致命错误 (readlink): 无效的 TITLE 行。id=%d, lid=%d。 id必须在[1,2]范围内, lid必须在[1,%d]范围内。\n", id + 1, lid + 1, BESTCLUSTERNUM);
                    printf("问题行: %s", line);
                    debug_flush();
                    fclose(pf);
                    exit(1);
                }
            } else {
                // TITLE行格式不匹配，可以忽略或警告
            }
        }
        else if(!strncmp(line,"REMARK  dock1",13)){
            if(hasgetfragname1==0){
                char *last_space = strrchr(line, ' ');
                if (last_space) {
                    char *path_start = last_space + 1;
                    
                    // *** 核心修改：在复制文件名之前，清理字符串 ***
                    trim_trailing_whitespace(path_start);

                    strncpy(fragname[0], path_start, 99);
                    fragname[0][99] = '\0'; // 确保安全
                } else {
                    fragname[0][0] = '\0'; // 未找到，设为空字符串
                }
                fragscore[0]=atof(line+16); // 假设分数位置固定
                hasgetfragname1=1;
            }
            if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                tertyp[id][lid][0]=atoi(line+13); // 假设tertyp位置固定
            }
        }
        else if(!strncmp(line,"REMARK  dock2",13)){
            if(hasgetfragname2==0){
                char *last_space = strrchr(line, ' ');
                if (last_space) {
                    char *path_start = last_space + 1;
                    
                    // *** 核心修改：在复制文件名之前，清理字符串 ***
                    trim_trailing_whitespace(path_start);

                    strncpy(fragname[1], path_start, 99);
                    fragname[1][99] = '\0'; // 确保安全
                } else {
                    fragname[1][0] = '\0'; // 未找到，设为空字符串
                }
                fragscore[1]=atof(line+16);
                hasgetfragname2=1;
            }
            if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                tertyp[id][lid][1]=atoi(line+13);
            }
        }
        else if(!strncmp(line,"REMARK  link",12)){
            if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                getfname(line,linkname[id][lid]);
            }
        }
        else if(!strncmp(line,"REMARK  score",13)){
            if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                linkscore[id][lid]=atof(line+14);
            }
        }
        if(!strncmp(line,"ATOM",4)){
            if(rn == 0 || strncmp(line+17, p[id][lid].res[rn-1].nam, 9) != 0){
                 if (rn >= MAX_RES_IN_FRAG) {
                    printf("致命错误 (readlink): 片段 id=%d lid=%d 的残基数过多 (%d)。最大允许 %d。\n", id, lid, rn, MAX_RES_IN_FRAG);
                    debug_flush();
                    fclose(pf);
                    exit(1);
                }
                if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                    strncpy(p[id][lid].res[rn].nam,line+17,9);
                    p[id][lid].res[rn].nam[9]='\0';
                    p[id][lid].res[rn].start=n;
                    p[id][lid].res[rn].sstart=n+4;
                    if(rn!=0) p[id][lid].res[rn-1].end=n-1;
                    rn++;
                } else {
                    printf("警告 (readlink): ATOM行在无效的 id/lid 区域被发现。id=%d, lid=%d。行: %s", id, lid, line);
                    debug_flush();
                    continue;
                }
            }
            
            if (n >= 100) {
                 printf("致命错误 (readlink): 片段 id=%d lid=%d 的原子数过多 (%d)。最大允许 100。\n", id, lid, n);
                 debug_flush();
                 fclose(pf);
                 exit(1);
            }
            
            if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                getxyz(line,p[id][lid].xyz[n]);
                strncpy(p[id][lid].nam[n],line+13,3);
                p[id][lid].nam[n][3]='\0';
                
                int current_res_index = rn-1;
                if (current_res_index < 0) {
                    printf("警告 (readlink): 在解析ATOMS时current_res_index为负值。id=%d, lid=%d, 行: %s", id, lid, line);
                    debug_flush();
                    continue; 
                }
                if(!strncmp(line+13,"N  ",3)){
                    p[id][lid].res[current_res_index].N=n;
                }
                else if(!strncmp(line+13,"CA ",3)){
                    p[id][lid].res[current_res_index].CA=n;
                }
                else if(!strncmp(line+13,"C  ",3)){
                    p[id][lid].res[current_res_index].C=n;
                }
                else if(!strncmp(line+13,"O  ",3)){
                    p[id][lid].res[current_res_index].O=n;
                }
                else if(!strncmp(line+13,"CB ",3)){
                    p[id][lid].res[current_res_index].sstart=n;
                }
                n++;
            }
        }
        if(!strncmp(line,"END",3)){
            if (id >= 0 && id < 2 && lid >= 0 && lid < BESTCLUSTERNUM) {
                p[id][lid].an=n;
                p[id][lid].resn=rn;
                if (rn > 0) p[id][lid].res[rn-1].end=n-1;
            } else {
                printf("警告 (readlink): END记录在无效的 id/lid 区域被发现。id=%d, lid=%d。行: %s", id, lid, line);
                debug_flush();
            }
            n=0;rn=0; 
        }
    }
    fclose(pf);
}

void cpx(float a[3],float b[3])
{
    int i;
    for(i=0;i<3;i++) a[i]=b[i];
}
void assemble(FRAGMENT *p1, FRAGMENT *l2, FRAGMENT *p2, FRAGMENT *l1, int t1, int t2, int tertyp1[2], int tertyp2[2],char chid,CYCPEP *cyc)
{
    int resn=0; // cyc->res 的残基索引
    int ir,ip;
    int sr,er; 
    int sa,ea; 
    int i;
    FRAGMENT *p;
    char resname[4]; 

    for(ip=0;ip<4;ip++){ 
        if(ip==0){ p=p1; sr=0;er=p1->resn;if(t1==1) sa=0; else sa=1; if(t2==1) ea=0; else ea=1;} 
        else  if(ip==1){ p=l2;if(t2==1){ sr=1;er=l2->resn-1;sa=0; ea=0;} else{sr=0;er=l2->resn;sa=2; ea=2;} } 
        else if(ip==2){ p=p2; sr=0;er=p2->resn;if(t2==1) sa=0; else sa=1; if(t1==1) ea=0; else ea=1;} 
        else if(ip==3){ p=l1; if(t1==1){ sr=1;er=l1->resn-1;sa=0; ea=0;} else{sr=0;er=l1->resn;sa=2; ea=2;}} 
        
        for(ir=sr;ir<er;ir++){ 
            if (resn >= MAX_RES_IN_PEP) {
                 printf("致命错误 (assemble): 尝试添加第 %d 个残基，但最大允许 %d。\n", resn, MAX_RES_IN_PEP);
                 debug_flush();
                 exit(1);
            }

            strncpy(resname, p->res[ir].nam, 3);
            resname[3] = '\0'; 

            if(!(sa==2&&ir==sr)){
            }
            if(ea==2&&ir==er-1){ 
                if(p==l2){ 
                    strncpy(resname, p2->res[0].nam, 3);
                    resname[3] = '\0';
                }
                else if(p==l1){ 
                    strncpy(resname, p1->res[0].nam, 3);
                    resname[3] = '\0';
                    resn=0; 
                }
            }
            sprintf(cyc->res[resn].resnam,"%s %c%4d",resname,chid,resn+1);
            cyc->res[resn].satn=0; 

            for(i=p->res[ir].start;i<=p->res[ir].end;i++){
                if(ir==sr){ 
                    if(sa!=0 && i==p->res[ir].N) continue; 
                    if(sa==2 && i==p->res[ir].CA) continue; 
                }
                if(ir==er-1){ 
                    if(ea==1){
                        if(i==p->res[ir].C||i==p->res[ir].O) continue;
                    }
                    if(ea==2){
                        if(i>=p->res[ir].CA) break; 
                    }
                }
                if(sa==2&&ir==sr){
                    if(i!=p->res[ir].C&&i!=p->res[ir].O) continue; 
                }
                
                if((p==p1)&&ir==sr&&tertyp1[0]==0&&i>=p->res[ir].sstart) break;
                if((p==p2)&&ir==sr&&tertyp2[0]==0&&i>=p->res[ir].sstart) break;
                if((p==p2)&&ir==er-1&&tertyp1[1]==0&&i>=p->res[ir].sstart) break;
                if((p==p1)&&ir==er-1&&tertyp2[1]==0&&i>=p->res[ir].sstart) break;

                if (cyc->res[resn].satn >= 10) { 
                     printf("致命错误 (assemble): 残基 %d ('%s') 的侧链原子过多 (%d)。最大允许 10。\n", resn, cyc->res[resn].resnam, cyc->res[resn].satn);
                     debug_flush();
                     exit(1);
                }

                if(!strcmp(p->nam[i],"N  ")) cpx(cyc->res[resn].N, p->xyz[i]);
                else if(!strcmp(p->nam[i],"CA ")) cpx(cyc->res[resn].CA, p->xyz[i]);
                else if(!strcmp(p->nam[i],"C  ")) cpx(cyc->res[resn].C, p->xyz[i]);
                else if(!strcmp(p->nam[i],"O  ")) cpx(cyc->res[resn].O, p->xyz[i]);
                else{ 
                    cpx(cyc->res[resn].side[cyc->res[resn].satn], p->xyz[i]);
                    strcpy(cyc->res[resn].nam[cyc->res[resn].satn], p->nam[i]); 
                    cyc->res[resn].nam[cyc->res[resn].satn][3]='\0'; 
                    cyc->res[resn].satn++;
                }
            }
            resn++;
            if(ir==er-1&&ea!=0){
                resn--;
            }
        }
    }
    cyc->resn = resn; 
}

void printpepstructure(FILE *pf,int id,char fragname[2][100], float dockscore[2], CYCPEP *cyc)
{
    int atno=1; 
    int ir,ia;
    int ih;

    fprintf(pf,"TITLE %2d\n",id+1);
    fprintf(pf,"REMARK  dock1 %s\n",fragname[0]);
    fprintf(pf,"REMARK  dock2 %s\n",fragname[1]);
    fprintf(pf,"REMARK  link1 %s %d\n",cyc->linkname[0], cyc->typ[0]);
    fprintf(pf,"REMARK  link2 %s %d\n",cyc->linkname[1], cyc->typ[1]);
    fprintf(pf,"REMARK  score %8.3f %8.3f %8.3f %8.3f %8.3f\n",dockscore[0]+dockscore[1]+cyc->linkscore[0]+cyc->linkscore[1],dockscore[0],dockscore[1],cyc->linkscore[0],cyc->linkscore[1]);
    fprintf(pf,"REMARK  resn %2d %3d %3d %3d %3d\n",cyc->resn,cyc->len[0],cyc->len[1],cyc->len[2],cyc->len[3]);
    fprintf(pf,"REMARK  backbone HB No. %2d\n",cyc->hbn);
    for(ih=0;ih<cyc->hbn;ih++){
        if (cyc->index[ih][0] < cyc->resn && cyc->index[ih][1] < cyc->resn) {
            fprintf(pf,"REMARK  HB %2d [%s]--[%s] %8.3f\n",ih+1,cyc->res[cyc->index[ih][0]].resnam,cyc->res[cyc->index[ih][1]].resnam, cyc->hbdis[ih]);
        } else {
            fprintf(pf,"WARNING: Invalid HB index detected: %d, %d for HB %d\n", cyc->index[ih][0], cyc->index[ih][1], ih+1);
            debug_flush();
        }
    }
    for(ir=0;ir<cyc->resn;ir++){
        fprintf(pf,"ATOM %6d  N   %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].N[0],cyc->res[ir].N[1],cyc->res[ir].N[2]); atno++;
        fprintf(pf,"ATOM %6d  CA  %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].CA[0],cyc->res[ir].CA[1],cyc->res[ir].CA[2]); atno++;
        fprintf(pf,"ATOM %6d  C   %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].C[0],cyc->res[ir].C[1],cyc->res[ir].C[2]); atno++;
        fprintf(pf,"ATOM %6d  O   %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].O[0],cyc->res[ir].O[1],cyc->res[ir].O[2]); atno++;
        for(ia=0;ia<cyc->res[ir].satn;ia++){
            fprintf(pf,"ATOM %6d %-4s%s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].nam[ia],cyc->res[ir].resnam,cyc->res[ir].side[ia][0],cyc->res[ir].side[ia][1],cyc->res[ir].side[ia][2]); atno++;
        }
    }

    fprintf(pf,"END\n");
}
int checkcollision(FRAGMENT *p1,FRAGMENT *p2)
{
        int ia,ja;
        float dis;

        for(ia=0;ia<p1->an;ia++){
                for(ja=0;ja<p2->an;ja++){
                        dis=distance2(p1->xyz[ia], p2->xyz[ja]);
                        if(dis<2.5*2.5){ 
                                return 0; // 发生碰撞
                        }
                }
        }
        return 1; // 无碰撞
}

#define MAXMODEL 25
void calHBond(CYCPEP *cyc)
{
    int pair[MAXHBN][2];
    int hbn,ih; 
    float hbdis[MAXHBN];
    int index[MAXHBN];

    addHN(cyc);
    hbn=dispair(cyc, pair, hbdis); 
    if(hbn > 0) {
      hbn=checkangle(cyc, hbn, pair, index);
    }
    cyc->hbn=hbn;
    for(ih=0;ih<hbn;ih++){
        cyc->index[ih][0]=pair[index[ih]][0];
        cyc->index[ih][1]=pair[index[ih]][1];
        cyc->hbdis[ih]=hbdis[index[ih]];
    }
}

int main(int argc, char *argv[])
{
    FRAGMENT (*p)[BESTCLUSTERNUM]; 
    FRAGMENT *p1, *p2;
    CYCPEP *cyc;

    char fragname[2][100];
    char linkname[2][BESTCLUSTERNUM][100];
    int typ[2][BESTCLUSTERNUM];
    float fragscore[2];
    float linkscore[2][BESTCLUSTERNUM];
    int tertyp[2][BESTCLUSTERNUM][2];
    int clusterN[2] = {0, 0};
    int il,jl;
    int col;
    int maxlinkres;
    int modn=0,imod; 
    FILE *pf;

    if (argc < 5) {
        printf("致命错误: 无效的参数数量。\n");
        printf("用法: %s <link_file> <output_pdb> <max_linker_res> <chain_id>\n", argv[0]);
        debug_flush();
        return 1;
    }
    
    p = malloc(2 * sizeof(*p)); 
    p1 = malloc(sizeof(FRAGMENT));
    p2 = malloc(sizeof(FRAGMENT));
    cyc = malloc(MAXMODEL * sizeof(CYCPEP));

    if (p == NULL || p1 == NULL || p2 == NULL || cyc == NULL) {
        printf("致命错误: 内存分配失败。\n");
        debug_flush();
        free(p);
        free(p1);
        free(p2);
        free(cyc);
        return 1;
    }
    
    fragname[0][0] = '\0';
    fragname[1][0] = '\0';

    maxlinkres=atoi(argv[3]);

    readlink(argv[1], p,typ,fragname,linkname,fragscore,linkscore,tertyp,clusterN);
    
    if(strlen(fragname[0]) == 0 || strlen(fragname[1]) == 0) {
        printf("\n致命错误: 未能从链接文件 '%s' 中正确读取片段文件名 (REMARK dock1/dock2)。\n\n", argv[1]);
        debug_flush();
        free(p); free(p1); free(p2); free(cyc);
        return 1;
    }

    readfrag(fragname[0], p1);
    readfrag(fragname[1], p2);

    for(il=0;il<clusterN[0];il++) {
        for(jl=0;jl<clusterN[1];jl++){
            if(modn>=MAXMODEL){
                printf("警告: 生成的模型数量 (%d) 已达到最大值 MAXMODEL (%d)，剩余组合将被忽略。\n", modn, MAXMODEL);
                debug_flush();
                goto end_loops; 
            }

            if(p[0][il].resn+p[1][jl].resn > maxlinkres) continue;
            
            col=checkcollision(&(p[1][jl]),&(p[0][il]));
            if(col==0) continue; 
            
            assemble(p1, &(p[1][jl]),p2,&(p[0][il]),typ[0][il],typ[1][jl],tertyp[0][il],tertyp[1][jl],argv[4][0],cyc+modn);
            
            strcpy(cyc[modn].linkname[0],linkname[0][il]); 
            strcpy(cyc[modn].linkname[1],linkname[1][jl]);
            cyc[modn].typ[0]=typ[0][il]; 
            cyc[modn].typ[1]=typ[1][jl];
            cyc[modn].linkscore[0]=linkscore[0][il]; 
            cyc[modn].linkscore[1]=linkscore[1][jl];
            cyc[modn].len[0]=p1->resn; 
            cyc[modn].len[1]=p[1][jl].resn;
            cyc[modn].len[2]=p2->resn; 
            cyc[modn].len[3]=p[0][il].resn;
            
            cyc[modn].resn=cyc[modn].len[0]+cyc[modn].len[1]+cyc[modn].len[2]+cyc[modn].len[3]-4; 
            
            modn++;
        }
    }
end_loops:; 

    if (modn == 0) {
        printf("\n---> 警告：程序未生成任何模型 <---\n");
        printf("这很可能是因为以下原因之一:\n");
        printf("1. 两个主要片段中，至少有一个没有在链接文件中找到对应的连接器(linker)簇。\n");
        printf("   根据链接文件 '%s' 的解析:\n", argv[1]);
        printf("   - 片段1 (dock1) 找到的连接器簇数量 (clusterN[0]): %d\n", clusterN[0]);
        printf("   - 片段2 (dock2) 找到的连接器簇数量 (clusterN[1]): %d\n", clusterN[1]);
        printf("   由于不满足此条件 (例如 %d * %d == 0)，没有模型被创建。\n", clusterN[0], clusterN[1]);
        printf("2. 所有可能的组合都因为碰撞检查或连接器残基数量限制而被跳过。\n");
        printf("   - <max_linker_res> 参数: %d\n", maxlinkres);
        printf("因此，输出文件 '%s' 将为空。\n\n", argv[2]);
        debug_flush();

        free(p);
        free(p1);
        free(p2);
        free(cyc);
        return 0; 
    }

    if((pf=fopen(argv[2],"w"))==NULL){
        printf("错误: 无法创建文件 %s\n",argv[2]);
        debug_flush();
        free(p); free(p1); free(p2); free(cyc);
        exit(1);
    }

    for(imod=0;imod<modn;imod++){
        calHBond(cyc+imod); 
        printpepstructure(pf, imod,fragname,fragscore, cyc+imod);
    }
    fclose(pf);

    free(p);
    free(p1);
    free(p2);
    free(cyc);

    printf("程序成功完成，共生成 %d 个模型。\n", modn);
    return 0; 
}