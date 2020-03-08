#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef struct
{
 unsigned char B,G,R;
}pixel;

typedef struct
{
  short T;
  long Size;
  short R1;
  short R2;
  long OffBits;
  long biSize;
  long biWidth;
  long biHeight;
  short biPlanes;
  short biBitCount;
  long biCompression;
  long biSizeImage;
  long biXPelsPerMeter;
  long biYPelsPerMeter;
  long biClrUsed;
  long biClrImportant;
} BMP_header;

typedef struct
{
 int poz;
 float corr;
 pixel cul;
}fereastra;

typedef struct
{
 int n;
 fereastra *f;
}ferestre;

typedef struct
{
 BMP_header header;
 pixel * p;
}imagine_BMP;

unsigned int xorshift32(unsigned int x)
{
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return x;
}

unsigned int *generare(int n , unsigned int R0)
{
 int i;
 unsigned int *v=malloc(sizeof(unsigned int)*n+1);
 v[0]=R0;
 for(i=1;i<=n;i++)
  v[i]=xorshift32(v[i-1]);

 return v;
}

imagine_BMP incarcare_BMP(char *cale)
{
 FILE *fin=fopen(cale,"rb");
 int i,n,j,height,width;

 imagine_BMP v;

 fread(&v.header.T, 2, 1, fin);
 fread(&v.header.Size, 4, 1, fin);
 fread(&v.header.R1, 2, 1, fin);
 fread(&v.header.R2, 2, 1, fin);
 fread(&v.header.OffBits, 4, 1, fin);
 fread(&v.header.biSize, 4, 1, fin);
 fread(&v.header.biWidth, 4, 1, fin);
 fread(&v.header.biHeight, 4, 1, fin);
 fread(&v.header.biPlanes, 2, 1, fin);
 fread(&v.header.biBitCount, 2, 1, fin);
 fread(&v.header.biCompression, 4, 1, fin);
 fread(&v.header.biSizeImage, 4, 1, fin);
 fread(&v.header.biXPelsPerMeter, 4, 1, fin);
 fread(&v.header.biYPelsPerMeter, 4, 1, fin);
 fread(&v.header.biClrUsed, 4, 1, fin);
 fread(&v.header.biClrImportant, 4, 1, fin);

 height=v.header.biHeight;
 width=v.header.biWidth;

 n=height*width+1;
 v.p=(pixel*)malloc(sizeof(pixel)*n+1);

 int padding;
    if(v.header.biWidth % 4 != 0)
        padding = 4 - (3 * v.header.biWidth) % 4;
    else
        padding = 0;

 for(i=height-1;i>=0;i--)
   {
    for(j=0;j<width;j++)
       fread(&v.p[(i*width)+j],sizeof(pixel),1,fin);

    fseek(fin,padding,SEEK_CUR);
   }

 fclose(fin);
 return v;
}

void salvare_BMP(char *cale , imagine_BMP v)
{
 FILE *fout=fopen(cale,"wb+");
 int i,j,k;

 fwrite(&v.header.T, 2, 1, fout);
 fwrite(&v.header.Size, 4, 1, fout);
 fwrite(&v.header.R1, 2, 1, fout);
 fwrite(&v.header.R2, 2, 1, fout);
 fwrite(&v.header.OffBits, 4, 1, fout);
 fwrite(&v.header.biSize, 4, 1, fout);
 fwrite(&v.header.biWidth, 4, 1, fout);
 fwrite(&v.header.biHeight, 4, 1, fout);
 fwrite(&v.header.biPlanes, 2, 1, fout);
 fwrite(&v.header.biBitCount, 2, 1, fout);
 fwrite(&v.header.biCompression, 4, 1, fout);
 fwrite(&v.header.biSizeImage, 4, 1, fout);
 fwrite(&v.header.biXPelsPerMeter, 4, 1, fout);
 fwrite(&v.header.biYPelsPerMeter, 4, 1, fout);
 fwrite(&v.header.biClrUsed, 4, 1, fout);
 fwrite(&v.header.biClrImportant, 4, 1, fout);

 int padding,p;
    if(v.header.biWidth % 4 != 0)
        padding = 4 - (3 * v.header.biWidth) % 4;
    else
        padding = 0;
 pixel pad;
 pad.B=pad.G=pad.R=0;

 for(i=v.header.biHeight-1;i>=0;i--)
  {for(j=0;j<v.header.biWidth;j++)
     fwrite(&v.p[(i*v.header.biWidth)+j],sizeof(pixel),1,fout);
   if(padding!=0)
   for(k=1;k<=p;k++)
     fwrite(&pad,sizeof(pixel),1,fout);
  }
 fclose(fout);
}

unsigned int *Durstenfeld(unsigned int *R , int n)
{
 unsigned int k,aux,r;
 unsigned int *p=malloc(sizeof(int)*n+1);
 for(k=0;k<n;k++)
     p[k]=k;

 for(k=n-1;k>=1;k--)
 {
  r=R[k]%(k+1);
  aux=p[r];
  p[r]=p[k];
  p[k]=aux;
 }

 return p;

}

pixel XOR_pixel_int(pixel p , unsigned int x)
{
 char *byte=&x;
 pixel y;

 y.B=p.B^byte[0];
 y.G=p.G^byte[1];
 y.R=p.R^byte[2];

 return y;
}

pixel *XOR_imagine(pixel *c , unsigned int *R , unsigned int SV , int n)
{
 int i;

 c[0]=XOR_pixel_int(XOR_pixel_int(c[0],SV),R[n]);

 for(i=1;i<n;i++)
    {
     c[i].B=c[i-1].B^c[i].B;
     c[i].G=c[i-1].B^c[i].G;
     c[i].R=c[i-1].B^c[i].R;

     c[i]=XOR_pixel_int(c[i],R[i+n]);
    }

 return c;
}

pixel *XOR_imagine_inv(pixel *c , unsigned int *R , unsigned int SV , int n)
{
 int i;


 for(i=n-1;i>1;i--)
    {
     c[i].B=c[i-1].B^c[i].B;
     c[i].G=c[i-1].B^c[i].G;
     c[i].R=c[i-1].B^c[i].R;

     c[i]=XOR_pixel_int(c[i],R[i+n]);
    }

 c[0]=XOR_pixel_int(XOR_pixel_int(c[0],SV),R[n]);

 return c;
}

pixel *permutare_pixeli(pixel *v, int *p, int n)
{
 pixel *v_perm=malloc(sizeof(pixel)*n+100);
 int i;

 for(i=0;i<n;i++)
  v_perm[i]=v[p[i]];
 return v_perm;
}

void criptare_BMP(char *imagine_initiala , char *imagine_criptata , char *secret_key)
{
 unsigned int i,n,*R,*p,*c,SV,R0;
 imagine_BMP v;

 FILE *secret=fopen(secret_key,"r");
 fscanf(secret,"%d %d" , &R0, &SV);

 v=incarcare_BMP(imagine_initiala);
 n=v.header.biHeight*v.header.biWidth;

 R=generare(n*2,R0);

 p=Durstenfeld(R,n);

 v.p=permutare_pixeli(v.p,p,n);

 v.p=XOR_imagine(v.p,R,SV,n);

 salvare_BMP(imagine_criptata,v);

 free(R);
 free(p);
 free(v.p);

 fclose(secret);
}

void decriptare_BMP(char *imagine_initiala , char *imagine_decriptata , char *secret_key)
{
 unsigned int i,n,*R,*p,*c,SV,R0;
 imagine_BMP v;

 FILE *secret=fopen(secret_key,"r");
 fscanf(secret,"%d %d" , &R0, &SV);

 v=incarcare_BMP(imagine_initiala);
 n=v.header.biHeight*v.header.biWidth;

 R=generare(n*2,R0);

 p=Durstenfeld(R,n);

 v.p=XOR_imagine_inv(v.p,R,SV,n);

 int *pi=malloc(sizeof(int)*n+1);

 for(i=0;i<n;i++)
   pi[p[i]]=i;

 v.p=permutare_pixeli(v.p,pi,n);

 salvare_BMP(imagine_decriptata,v);

 free(R);
 free(p);
 free(pi);
 free(v.p);
 fclose(secret);
}

void chi_patrat(char *cale)
{
 imagine_BMP v;

 v=incarcare_BMP(cale);

 int n=v.header.biHeight*v.header.biWidth;

 int *frec_R=malloc(sizeof(int)*257);
 int *frec_G=malloc(sizeof(int)*257);
 int *frec_B=malloc(sizeof(int)*257);

 int i;

  for(i=0;i<=255;i++)
   {
    frec_R[i]=0;
    frec_G[i]=0;
    frec_B[i]=0;
   }

 for(i=0;i<n;i++)
   {
    frec_R[v.p[i].R]++;
    frec_G[v.p[i].G]++;
    frec_B[v.p[i].B]++;
   }

 float frec_med=n/256;


 float sum_R=0,sum_B=0,sum_G=0;

 for(i=0;i<=255;i++)
 {
  sum_R=sum_R+(((frec_R[i]-frec_med)*(frec_R[i]-frec_med))/frec_med);
  sum_G=sum_G+(((frec_G[i]-frec_med)*(frec_G[i]-frec_med))/frec_med);
  sum_B=sum_B+(((frec_B[i]-frec_med)*(frec_B[i]-frec_med))/frec_med);
 }

 printf("R:%f\nG:%f\nB:%f\n" , sum_R, sum_G, sum_B);

 free(frec_R);
 free(frec_G);
 free(frec_B);
}

void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
   FILE *fin, *fout;
   unsigned int dim_img, latime_img, inaltime_img;
   unsigned char pRGB[3], header[54], aux;

   printf("nume_fisier_sursa = %s \n",nume_fisier_sursa);

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

   fout = fopen(nume_fisier_destinatie, "wb+");

   fseek(fin, 2, SEEK_SET);
   fread(&dim_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);

   //copiaza octet cu octet imaginea initiala in cea noua
	fseek(fin,0,SEEK_SET);
	unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);

	//calculam padding-ul pentru o linie
	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    printf("padding = %d \n",padding);

	fseek(fout, 54, SEEK_SET);
	int i,j;
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{
			//citesc culorile pixelului
			fread(pRGB, 3, 1, fout);
			//fac conversia in pixel gri
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout,padding,SEEK_CUR);
	}
	fclose(fout);
}

float intensitate_medie(pixel *v , int n)
{
 int i,sum=0;
 for(i=0;i<n;i++)
    sum+=v[i].R;
 return sum/n;
}

float deviatie(pixel *v, int n , float med)//med este transmis ca parametru pentru a creste eficienta
{
 int i;
 float sum=0;

 for(i=0;i<n;i++)
    sum+=((v[i].R-med)*(v[i].R-med));
 sum=sum/(n-1);
 return sqrt(sum);

}

float corr(pixel *S , pixel *f ,float S_med ,int n)//S_med este transmis ca parametru pentru a creste eficienta
{
 int i;
 float sum=0;
 float f_med=intensitate_medie(f,n);
 float S_dev=(deviatie(S,n,S_med));
 for(i=0;i<n;i++)
  sum+=((f[i].R-f_med)*(S[i].R-S_med)/(S_dev*deviatie(f,n,f_med)));

 return sum/n;
}

ferestre template_matching(imagine_BMP v, imagine_BMP S, float pS, pixel culoare)
{
 int i,j,k,sol=0;
 int ns=S.header.biHeight*S.header.biWidth;
 int nv=v.header.biHeight*v.header.biWidth;
 float S_med=intensitate_medie(S.p,ns);
 imagine_BMP f_cur;
 f_cur.p=(pixel*)malloc(sizeof(pixel)*ns+1);//fereastra curenta
 ferestre fer;
 fer.f=(ferestre *)malloc(sizeof(fereastra)*nv);

 //parcurgem vectorul liniarizat ca si cum ar fi o matrice
 for(i=0;i<nv-(15*v.header.biWidth);i++)//evitam cazul in care sablonus iese din imagine in jos
 {
  if((v.header.biWidth-(i%v.header.biWidth))>S.header.biWidth)//evitam cazul in care sablon iese din imagine in dreapta
  {
   for(j=0;j<S.header.biHeight;j++)
   {
    for(k=0;k<S.header.biWidth;k++)
    {
     int x1=k+j*S.header.biWidth;
     int x2=(i+k)+(j*v.header.biWidth);
     f_cur.p[x1]=v.p[x2];
    }
   }

   float x=corr(S.p,f_cur.p,S_med,ns);//calculam corelatia
   if(x>pS)//daca corelatia este peste un prag, retinem fereastra
   {
    fer.f[sol].poz=i;
    fer.f[sol].cul=culoare;
    fer.f[sol].corr=x;
    sol++;
   }
  }
 }
 fer.f=realloc(fer.f,sizeof(fereastra)*(sol+1));//realocam fereastra pentru a ocupa mai putina memorie
 fer.n=sol;

 free(f_cur.p);
 return fer;
}

ferestre template_matching_2(imagine_BMP v , char ** cale)
{
 int i,j;
 ferestre fer,det;
 imagine_BMP S;
 det.f=(ferestre*)malloc(sizeof(fereastra)*v.header.biHeight*v.header.biWidth);
 det.n=0;

 pixel *culoare=malloc(sizeof(pixel)*11);
 culoare[0].R=255; culoare[0].G=0; culoare[0].B=0;
 culoare[1].R=255; culoare[1].G=255; culoare[1].B=0;
 culoare[2].R=0; culoare[2].G=255; culoare[2].B=0;
 culoare[3].R=0; culoare[3].G=255; culoare[3].B=255;
 culoare[4].R=255; culoare[4].G=0; culoare[4].B=255;
 culoare[5].R=0; culoare[5].G=0; culoare[5].B=255;
 culoare[6].R=192; culoare[6].G=192; culoare[6].B=192;
 culoare[7].R=255; culoare[7].G=140; culoare[7].B=0;
 culoare[8].R=128; culoare[8].G=0; culoare[8].B=128;
 culoare[9].R=120; culoare[9].G=0; culoare[9].B=0;

 for(i=0;i<10;i++)
 {
  S=incarcare_BMP(cale[i]);
  fer=template_matching(v,S,0.5,culoare[i]);

  for(j=0;j<fer.n;j++)
     det.f[j+det.n]=fer.f[j];
  det.n+=fer.n;

  free(fer.f);
 }

 det.f=realloc(det.f,sizeof(fereastra)*det.n+1);

 return det;
}

imagine_BMP colorare_contur(imagine_BMP v , imagine_BMP S, int f , pixel culoare)
{
 int i;
 for(i=0;i<S.header.biWidth;i++)
  {
   v.p[f+i]=culoare;
   v.p[f+i+(S.header.biHeight*v.header.biWidth)]=culoare;
  }
 for(i=0;i<S.header.biHeight;i++)
  {
   v.p[f+i*v.header.biWidth]=culoare;
   v.p[f+i*v.header.biWidth+S.header.biWidth]=culoare;
  }

return v;
}

int comp_ferestre(const void *a ,const void*b)
{
 fereastra x= *(fereastra*) a;
 fereastra y= *(fereastra*) b;
 if(x.corr>y.corr) return -1;
 if(x.corr<y.corr) return 1;
return 0;
}

float suprapunere(int f1 , int f2 , BMP_header header_S , BMP_header header_i)
{

 int x1,x2,y1,y2,i,arie_s,arie_inter,aux;

 if(f1>f2)//eliminam 2 cazuri pentru f1>f2
 {
  aux=f1;
  f1=f2;
  f2=aux;
 }

 if((f1%500)<=(f2%500))//cazul in care f1 e mai la stanga de f2
 {
  x1=(f1+header_S.biHeight*header_i.biWidth-1)/header_i.biWidth;
  x2=f2/header_i.biWidth;
  y1=(f1+header_S.biWidth-1)%header_i.biWidth;
  y2=f2%header_i.biWidth;
  arie_inter=(x1-x2+1)*(y1-y2+1);
  arie_s=header_S.biHeight*header_S.biWidth*2;
 }

 if((f1%500)>(f2%500))//cazul in care f1 e mai la dreapta de f2
 {
  x1=f1/header_i.biWidth;
  x2=(f2+header_S.biHeight*header_i.biWidth-1)/header_i.biWidth;
  y1=f1%header_i.biWidth;
  y2=(f2+header_S.biWidth-1)%header_i.biWidth;
  arie_inter=(x1-x2+1)*(y1-y2+1);
  arie_s=header_S.biHeight*header_S.biWidth*2;
 }


 float x=(float)arie_inter/(arie_s-arie_inter);
 return x;
}

imagine_BMP eliminarea_non_maximelor(imagine_BMP I , ferestre det , char * cale_S)
{

 int i,j;
 imagine_BMP S=incarcare_BMP(cale_S);
 qsort(det.f,det.n,sizeof(fereastra),comp_ferestre);

 for(i=0;i<det.n-1;i++)
 {
  for(j=i+1;j<det.n;j++)
     if(abs(det.f[i].poz-det.f[j].poz)<S.header.biHeight*I.header.biWidth)//cazul in care nu au nici o zona suprapusa
     {
      if(suprapunere(det.f[i].poz,det.f[j].poz,S.header,I.header)>0.2)
         det.f[j].corr=0;
     }
 }

 for(i=0;i<det.n;i++)
   if(det.f[i].corr!=0)
    I=colorare_contur(I,S,det.f[i].poz,det.f[i].cul);

 return I;

}

int main()
{
 imagine_BMP I,I_g;//imaginea si imaginea grayscale
 int i,j;
 ferestre det;//structura in care retine detectiile

 char *cale_imagine,*cale_imagine_criptata,*cale_imagine_decriptata,*cale_cheie_secreta,*cale_imagine_numere,*cale_imagine_numere_g,*cale_imagine_numere_desenata;
 char **cale_sablon=malloc(sizeof(char*)*11);
 char **cale_sablon_g=malloc(sizeof(char*)*11);
 imagine_BMP *sablon=malloc(sizeof(imagine_BMP)*11);
 cale_imagine=(char*) malloc(sizeof(char)*100);
 cale_imagine_criptata=(char*) malloc(sizeof(char)*100);
 cale_imagine_decriptata=(char*) malloc(sizeof(char)*100);
 cale_cheie_secreta=(char*) malloc(sizeof(char)*100);
 cale_imagine_numere=(char*) malloc(sizeof(char)*100);
 cale_imagine_numere_g=(char*) malloc(sizeof(char)*100);
 cale_imagine_numere_desenata=(char*) malloc(sizeof(char)*100);

 for(i=0;i<10;i++)
   {
    cale_sablon[i]=(char*) malloc(sizeof(char*)*100);//cale sabloane
    cale_sablon_g[i]=(char*) malloc(sizeof(char*)*100);//cale sabloane grayscale
   }

 cale_sablon_g[0]="0_g.bmp";
 cale_sablon_g[1]="1_g.bmp";
 cale_sablon_g[2]="2_g.bmp";
 cale_sablon_g[3]="3_g.bmp";
 cale_sablon_g[4]="4_g.bmp";
 cale_sablon_g[5]="5_g.bmp";
 cale_sablon_g[6]="6_g.bmp";
 cale_sablon_g[7]="7_g.bmp";
 cale_sablon_g[8]="8_g.bmp";
 cale_sablon_g[9]="9_g.bmp";
 cale_imagine_numere_g="imagine_numere_g.bmp";


//partea intai a proiectului

 printf("Cale imagine:");
 scanf("%s" , cale_imagine);

 printf("Cale imagine criptata:");
 scanf("%s" , cale_imagine_criptata);

 printf("Cale imagine decriptata:");
 scanf("%s" , cale_imagine_decriptata);

 printf("Cale cheie secreta:");
 scanf("%s" , cale_cheie_secreta);
 printf("\n");

 criptare_BMP(cale_imagine,cale_imagine_criptata,cale_cheie_secreta);
 decriptare_BMP(cale_imagine_criptata,cale_imagine_decriptata,cale_cheie_secreta);

 printf("Rezultatele testului chi-patrat pentru imaginea initiala:\n");
 chi_patrat(cale_imagine);
 printf("\n");
 printf("Rezultatele testului chi-patrat pentru imaginea criptata:\n");
 chi_patrat(cale_imagine_criptata);
 printf("\n");


 //partea a doua a proiectului

 printf("Partea a doua a proiectului:\n");

 printf("Cale imagine numere:");
 scanf("%s" , cale_imagine_numere);
 printf("Cale imagine numere desenata:");
 scanf("%s" , cale_imagine_numere_desenata);

 grayscale_image(cale_imagine_numere,cale_imagine_numere_g);
 I=incarcare_BMP(cale_imagine_numere);
 I_g=incarcare_BMP(cale_imagine_numere_g);

 for(i=0;i<10;i++)
 {
  printf("Cale sablon cifra %d: " , i);
  scanf("%s" , cale_sablon[i]);
 }

 for(i=0;i<10;i++)
  grayscale_image(cale_sablon[i],cale_sablon_g[i]);

 det=template_matching_2(I_g,cale_sablon_g);

 I=eliminarea_non_maximelor(I,det,cale_sablon_g[0]);

 salvare_BMP(cale_imagine_numere_desenata , I );



 free(I.p);
 free(det.f);
 free(cale_imagine);
 free(cale_imagine_criptata);
 free(cale_imagine_decriptata);
 free(cale_imagine_numere);
 free(cale_imagine_numere_desenata);
 free(cale_imagine_numere_g);
 free(cale_cheie_secreta);
 free(cale_sablon);
 free(cale_sablon_g);


    return 0;
}
