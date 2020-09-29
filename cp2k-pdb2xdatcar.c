/*BY FILIPE MATUSALEM, AUG 2020     filipematus@gmail.com */
/*Program to convert CP2K PDB trajectory file to VASP XDATCAR format*/
/*Compilation: g++ -o cp2k-pdb2xdatcar cp2k-pdb2xdatcar.c*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846


int main(int argc, char *argv[])
{
FILE *pdb,*xyz,*xdatcar;
float vecx[3],vecy[3],vecz[3],alpha,beta,gamma,x,y,z,a,b,c;
int i,j,k,natoms,nspecies,nsteps,ntype[10];
char str[150],str1[150],ch,species[10][10];


 if( argc < 2 ){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");

printf("Enter the name of the cp2k trajectory pdb file \n\n");

 exit(0);}


printf("-----------------------------------------------------------------------------------------------\n\n");



strcpy(str1,argv[1]);

pdb = fopen(str1,"r"); /* Arquivo ASCII, para leitura */
if(!pdb)
{
printf( "Error opening argument 1 file\n");
exit(0);
}

xdatcar = fopen("XDATCAR","w"); /* Arquivo ASCII, para escrita */
if(!xdatcar)
{
printf( "Error creating xdatcar file\n");
exit(0);
}


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

natoms=0;
do
{fscanf(pdb,"%s",str1);  if(strcmp(str1,"ATOM")==0)natoms++;  }                                   
while(strcmp(str1,"END")!=0);

printf("No atoms %d\n",natoms);
rewind(pdb);


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",species[0]);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

j=k=1;
for(i=0;i<natoms-1;i++){
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",species[j]);
k++;
if(strcmp(species[j-1],species[j])!=0){ntype[j-1]=k-1;k=1;j++;}

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');
}
ntype[j-1]=k;
nspecies=j;

printf("No. Species %d \n\n",nspecies);
printf("Specie  Number\n");
for(i=0;i<nspecies;i++){
printf("   %s      %d \n",species[i],ntype[i]);
}
rewind(pdb);

nsteps=0;
while (fscanf(pdb,"%s",str1) != EOF){            /*conta steps*/
if(strcmp(str1,"CRYST1")==0)nsteps++;                      
}

printf("No. steps = %d\n",nsteps);
rewind(pdb);

float M[3][3],V;


for(i=0;i<nsteps;i++){
do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

fscanf(pdb,"%f",&a);
fscanf(pdb,"%f",&b);
fscanf(pdb,"%f",&c);
fscanf(pdb,"%f",&alpha);
fscanf(pdb,"%f",&beta);
fscanf(pdb,"%f",&gamma);


alpha=PI*alpha/180;
beta=PI*beta/180;
gamma=PI*gamma/180;

//conversion from fractional coordinates system to cartesian
//http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm

V=a*b*c*sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)+2*cos(alpha)*cos(beta)*cos(gamma));

//transformation matrix M

M[0][0]=a;
M[0][1]=b*cos(gamma);
M[0][2]=c*cos(beta);

M[1][0]=0;
M[1][1]=b*sin(gamma);
M[1][2]=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);

M[2][0]=0;
M[2][1]=0;
M[2][2]=V/(a*b*sin(gamma));

//Vector_x = M * (1 0 0)^t
//Vector_y = M * (0 1 0)^t
//Vector_z = M * (0 0 1)^t

for(j=0;j<3;j++){
vecx[j]=M[j][0]*1+M[j][1]*0+M[j][2]*0;
vecy[j]=M[j][0]*0+M[j][1]*1+M[j][2]*0;
vecz[j]=M[j][0]*0+M[j][1]*0+M[j][2]*1;}


fprintf(xdatcar,"Converted from cp2k PDB\n");
fprintf(xdatcar,"           1\n");
fprintf(xdatcar,"     %10.6f %10.6f %10.6f\n",vecx[0],vecx[1],vecx[2]);
fprintf(xdatcar,"     %10.6f %10.6f %10.6f\n",vecy[0],vecy[1],vecy[2]);
fprintf(xdatcar,"     %10.6f %10.6f %10.6f\n",vecz[0],vecz[1],vecz[2]);

for(j=0;j<nspecies;j++){
fprintf(xdatcar,"   %s",species[j]);
}
fprintf(xdatcar,"\n");

for(j=0;j<nspecies;j++){
fprintf(xdatcar,"   %d",ntype[j]);
}
fprintf(xdatcar,"\n");

fprintf(xdatcar,"Cartesian configuration=     %d\n",i+1);
for(j=0;j<natoms;j++){
do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

fscanf(pdb,"%s",str1);fscanf(pdb,"%s",str1);fscanf(pdb,"%s",str1);
fscanf(pdb,"%f",&a);fprintf(xdatcar,"   %10.6f",a);
fscanf(pdb,"%f",&a);fprintf(xdatcar,"   %10.6f",a);
fscanf(pdb,"%f",&a);fprintf(xdatcar,"   %10.6f\n",a);
}






}


printf("\n");
printf("XDATCAR written!! \n");
printf("-----------------------------------------------------------------------------------------------\n\n");



fclose(pdb);
fclose(xdatcar);
}
