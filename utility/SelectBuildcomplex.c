#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXMODEL 25

// selectmod 函数保持不变
int selectmod(char *fn, float dG, int hbn, float tot, int bestN, int *index)
{
    FILE *pf;
    float score[MAXMODEL];
    int HB[MAXMODEL];
    int modn=0,imod;
    char line[100];
    int i,j,k;
    int maxi;

    if((pf=fopen(fn,"r"))==NULL){
        printf("ERROR: Can not open protein structure file1 %s\n",fn);
        exit(0);
    }
    while(fgets(line, 100,pf)){
        if(!strncmp(line,"TITLE",5)){
            if(modn < MAXMODEL) {
                score[modn]=0;
                HB[modn]=0;
                modn++;
            }
        }
        else if(!strncmp(line,"REMARK",6)){
            if (modn > 0) {
                if(!strncmp(line+8,"score",5)){
                    score[modn-1]=atof(line+13);
                }
                else if(!strncmp(line+8,"backbone HB No",14)){
                    HB[modn-1]=atoi(line+23);
                }
            }
        }
    }
    fclose(pf);

    i=0;
    for(imod=0;imod<modn;imod++){
        if(score[imod]<=dG&&HB[imod]>=hbn&&(score[imod]-HB[imod]*1.5)<=tot){
            if(i<bestN){
                index[i]=imod;
                i++;
            }
            else{
                index[bestN]=imod; maxi=bestN;
                for(j=0;j<bestN;j++){
                    if((score[maxi]-HB[maxi]*1.5)<(score[j]-HB[j]*1.5)) maxi=j;
                }
                if(maxi!=bestN) index[maxi]=imod;
            }
        }
    }
    return i;
}

// writetarget 函数保持不变
void writetarget(FILE *f1, char *target)
{
    FILE *pf;
    char line[100];
    
    if((pf=fopen(target,"r"))==NULL){
        printf("ERROR: Can not open protein structure file2 %s\n",target);
        exit(0);
    }
    while(fgets(line, 100,pf)){
        if(!strncmp(line,"ATOM",4)){
            line[56]='\n';
            fprintf(f1,"%s",line);
        }
    }
    fclose(pf);
}

// selected 函数保持不变
int selected(int modn, int N, int *index)
{
    int i;
    for(i=0;i<N;i++) if(modn==index[i]) return 1;
    return 0;
}


// --- MODIFIED: getname 函数被大大简化 ---
// 它现在接收一个 `prefix` 参数，而不是 `fn`，并直接用它来构建文件名。
void getname(char *cmplx, char *prefix, char *dir, int modn)
{
    // 之前从 fn 提取基础名的复杂逻辑被移除。
    sprintf(cmplx, "%s/%s_%d.pdb", dir, prefix, modn + 1);
}

// --- MODIFIED: buildcomplex 函数签名增加了一个 `prefix` 参数 ---
void buildcomplex(char *fn, int N, int *index, char *dir, char *target, char *prefix)
{
    FILE *pf,*cf;

    int modn=0;
    char line[100];
    char cmplx[100];
    int start=0;
    int atmstart=0;

    if((pf=fopen(fn,"r"))==NULL){
        printf("ERROR: Can not open protein structure file3 %s\n",fn);
        exit(0);
    }
    while(fgets(line, 100,pf)){
        if(!strncmp(line,"TITLE",5)){
            if(selected(modn,N,index)==1){
                // --- MODIFIED: 调用 getname 时传入 prefix ---
                getname(cmplx, prefix, dir, modn);
                if((cf=fopen(cmplx,"w"))==NULL){
                    // --- MODIFIED: 修正了错误信息中的变量名 ---
                    printf("ERROR: Can not create protein structure file4 %s\n",cmplx);
                    exit(0);
                }
                start=1;
            }
            modn++;
        }
        else if(start){
            if(!strncmp(line,"ATOM",4)){
                if(atmstart==0){
                    writetarget(cf,target);
                    fprintf(cf,"TER\n");
                    atmstart=1;
                }   
            }
            fprintf(cf,"%s",line);
            if(!strncmp(line,"END",3)){
                fclose(cf);atmstart=0;start=0;
            }
        }
    }
    fclose(pf);
}

// --- MODIFIED: main 函数现在处理 8 个参数 ---
int main(int argc, char *argv[])
{
    // --- MODIFIED: 参数数量检查变为 9 (程序名 + 8个参数) ---
    if (argc != 9) {
        // --- MODIFIED: 更新了用法说明 ---
        printf("Usage: %s <input_pdb> <dG_cutoff> <min_hbn> <tot_cutoff> <best_n> <output_dir> <target_pdb> <output_prefix>\n", argv[0]);
        return 1;
    }

    int index[MAXMODEL];
    float dG=atof(argv[2]);
    int hbn=atoi(argv[3]);
    float tot=atof(argv[4]);
    int bestN=atoi(argv[5]);
    int N;

    N=selectmod(argv[1], dG, hbn, tot, bestN, index);

    // --- MODIFIED: 调用 buildcomplex 时传入新的第8个参数 ---
    buildcomplex(argv[1], N, index, argv[6], argv[7], argv[8]);
    
    printf("====> %2d complex model built for %s with prefix %s\n", N, argv[1], argv[8]);

    return 0;
}