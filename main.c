#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct pixel{
    unsigned char R;
    unsigned char G;
    unsigned char B;
}pixel;

typedef struct corelatie{
    unsigned int linie;
    unsigned int col;
    double corr;
    unsigned int cifra;
}corelatie;

unsigned int * xorshift32(unsigned int seed,unsigned int n)
{

    unsigned int r,k;
    unsigned int *vec;
    vec=(unsigned int*)malloc(n*sizeof(unsigned int));
    r=seed;
    for(k=1;k<n;k++)
    {
        r=r^r<<13;
        r=r^r>>17;
        r=r^r<<5;
        vec[k]=r;
    }
    return vec;

}

void citire_imagine(char*nume_sursa,pixel **L,unsigned char **header)
{
    FILE*fin;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char pRGB[3];
    *header=(unsigned char*)malloc(54*sizeof(unsigned char));
    //printf("nume_fisier_sursa = %s \n",nume_sursa);

    fin = fopen(nume_sursa, "rb");
    if(fin == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }
    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    *L=(pixel*)malloc((inaltime_img*latime_img)*sizeof(pixel));
    int padding;
    if(latime_img % 4 != 0)
        padding =4- (3 * latime_img) % 4;
    else
        padding = 0;

    //printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);
    //printf("Padding = %d \n",padding);

    fseek(fin,0,SEEK_SET);
    unsigned int i,j;
    for(i=0;i<54;i++)
        fread((*header+i),sizeof(unsigned char),1,fin);
    fseek(fin, 54, SEEK_SET);
    for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {

            fread(pRGB, 3, 1, fin);
            (*L+((inaltime_img-1-i)*latime_img+j))->B=pRGB[0];
            (*L+((inaltime_img-1-i)*latime_img+j))->G=pRGB[1];
            (*L+((inaltime_img-1-i)*latime_img+j))->R=pRGB[2];

        }


        fseek(fin,padding,SEEK_CUR);

    }

    fclose(fin);
}
void afisare(char*nume_sursa,char*nume_destinatie,pixel *L,unsigned char *header)
{
    FILE*fin;
    FILE*fout;


    fin = fopen(nume_sursa, "rb");
    if(fin == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }


    fout = fopen(nume_destinatie, "wb");
    if(fout == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }
    unsigned int dim_img, latime_img, inaltime_img;


    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    int padding;
    if(latime_img % 4 != 0)
        padding =4-(3 * latime_img) % 4;
    else
        padding = 0;
    //printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);
    //printf("Padding = %d \n",padding);

    unsigned int i;
    fwrite(header,sizeof(unsigned char),54,fout);
    fflush(fout);

    unsigned int j;
    fseek(fin,54,SEEK_SET);
    for(i=0;i<inaltime_img;i+=1)
    {
        for(j=0;j<latime_img;j++)
        {

            fwrite(&L[(inaltime_img-1-i)*latime_img+j].B,1,1,fout);
            fwrite(&L[(inaltime_img-1-i)*latime_img+j].G,1,1,fout);
            fwrite(&L[(inaltime_img-1-i)*latime_img+j].R,1,1,fout);
            fflush(fout);
        }
        unsigned int k=0;
        while(k<padding)
        {
            unsigned char c='0';
            fwrite(&c,1,1,fout);
            k++;
        }
    }

    fclose(fin);
    fclose(fout);
}

pixel* criptare(char*nume_sursa,char*nume_destinatie,char*cheie_secreta)
{
    FILE*fin;
    unsigned int latime_img, inaltime_img,i;
    unsigned char*header;
    //printf("nume_fisier_sursa = %s \n",nume_sursa);

    fin = fopen(nume_sursa, "rb");
    if(fin == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    fclose(fin);

    unsigned int *R=(unsigned int*)malloc((2*inaltime_img*latime_img)*sizeof(unsigned int));
    unsigned int Ro,SV;
    fin=fopen(cheie_secreta,"r");
    fscanf(fin,"%u",&Ro);
    fscanf(fin,"%u",&SV);
    fclose(fin);
    R=xorshift32(Ro,2*inaltime_img*latime_img);
    R[0]=Ro;
    unsigned int *permutare;
    permutare=(unsigned int*)malloc((inaltime_img*latime_img)*sizeof(unsigned int));
    for(i=0;i<inaltime_img*latime_img;i++)
        permutare[i]=i;
    unsigned int r,cnt_R=1,aux;
    for(i=inaltime_img*latime_img-1;i>0;i--)
    {
        r=R[cnt_R]%(i+1);
        aux=permutare[r];
        permutare[r]=permutare[i];
        permutare[i]=aux;
        cnt_R++;
    }
    pixel*L,*L_permut;
    L=(pixel*)malloc((inaltime_img*latime_img)*sizeof(pixel));
    L_permut=(pixel*)malloc((inaltime_img*latime_img)*sizeof(pixel));

    citire_imagine(nume_sursa,&L,&header);
    free(header);
    for(i=0;i<inaltime_img*latime_img;i++)
    {

        L_permut[permutare[i]]=L[i];

    }

    pixel*C;
    C=(pixel*)malloc((inaltime_img*latime_img)*sizeof(pixel));
    unsigned char*pointer_octet,*pointer_R;
    pointer_octet=&SV;
    pointer_R=&R[inaltime_img*latime_img];
    C[0].B=L_permut[0].B^*(pointer_octet)^*(pointer_R);
    C[0].G=L_permut[0].G^*(pointer_octet+1)^*(pointer_R+1);
    C[0].R=L_permut[0].R^*(pointer_octet+2)^*(pointer_R+2);
    for(i=1;i<inaltime_img*latime_img;i++)
    {
        pointer_R=&R[inaltime_img*latime_img+i];
        C[i].B=L_permut[i].B^C[i-1].B^*(pointer_R);
        C[i].G=L_permut[i].G^C[i-1].G^*(pointer_R+1);
        C[i].R=L_permut[i].R^C[i-1].R^*(pointer_R+2);
    }

    free(L);
    free(L_permut);
    free(permutare);
    free(R);
    return C;

}

pixel* decriptare(char*nume_sursa,char*nume_destinatie,char*cheie_secreta)
{
    FILE*fin;
    unsigned int latime_img, inaltime_img,i;
    unsigned char*header;
    //printf("nume_fisier_sursa = %s \n",nume_sursa);

    fin = fopen(nume_sursa, "rb");
    if(fin == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    fclose(fin);

    unsigned int *R=(unsigned int*)malloc((2*inaltime_img*latime_img)*sizeof(unsigned int));
    unsigned int Ro,SV;
    fin=fopen(cheie_secreta,"r");
    fscanf(fin,"%u",&Ro);
    fscanf(fin,"%u",&SV);
    fclose(fin);
    R=xorshift32(Ro,2*inaltime_img*latime_img);
    R[0]=Ro;
    unsigned int *permutare,*permutare_inversa;
    permutare=(unsigned int*)malloc((inaltime_img*latime_img)*sizeof(unsigned int));
    permutare_inversa=(unsigned int*)malloc((inaltime_img*latime_img)*sizeof(unsigned int));
    for(i=0;i<inaltime_img*latime_img;i++)
        permutare[i]=i;
    unsigned int r,cnt_R=1,aux;
    for(i=inaltime_img*latime_img-1;i>0;i--)
    {
        r=R[cnt_R]%(i+1);
        aux=permutare[r];
        permutare[r]=permutare[i];
        permutare[i]=aux;
        cnt_R++;
    }
    for(i=0;i<inaltime_img*latime_img;i++)
    {
        permutare_inversa[permutare[i]]=i;
    }
    pixel*C,*C_permut;
    C=(pixel*)malloc((inaltime_img*latime_img)*sizeof(pixel));
    C_permut=(pixel*)malloc((inaltime_img*latime_img)*sizeof(pixel));

    citire_imagine(nume_sursa,&C,&header);
    free(header);

    pixel*D;
    D=(pixel*)malloc((inaltime_img*latime_img)*sizeof(pixel));
    unsigned char*pointer_octet,*pointer_R;
    pointer_octet=&SV;
    pointer_R=&R[inaltime_img*latime_img];
    C_permut[0].B=*(pointer_octet)^C[0].B^*(pointer_R);
    C_permut[0].G=*(pointer_octet+1)^C[0].G^*(pointer_R+1);
    C_permut[0].R=*(pointer_octet+2)^C[0].R^*(pointer_R+2);
    for(i=1;i<inaltime_img*latime_img;i++)
    {
        pointer_R=&R[inaltime_img*latime_img+i];
        C_permut[i].B=C[i-1].B^C[i].B^*(pointer_R);
        C_permut[i].G=C[i-1].G^C[i].G^*(pointer_R+1);
        C_permut[i].R=C[i-1].R^C[i].R^*(pointer_R+2);
    }
    for(i=0;i<inaltime_img*latime_img;i++)
    {
        D[permutare_inversa[i]]=C_permut[i];
    }

    free(C);
    free(C_permut);
    free(permutare);
    free(permutare_inversa);
    free(R);
    return D;

}
void chi(char*nume_sursa,char*cheie_secreta)
{
    unsigned int*frecventaB=(unsigned int*)calloc(256,sizeof(unsigned int));
    unsigned int*frecventaG=(unsigned int*)calloc(256,sizeof(unsigned int));
    unsigned int*frecventaR=(unsigned int*)calloc(256,sizeof(unsigned int));
    FILE*fin;
    unsigned int latime_img, inaltime_img,i;
    unsigned char *header;
    //printf("nume_fisier_sursa = %s \n",nume_sursa);

    fin = fopen(nume_sursa, "rb");
    if(fin == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    fclose(fin);
    pixel*liniarizare;
    citire_imagine(nume_sursa,&liniarizare,&header);
    free(header);

    for(i=0;i<latime_img*inaltime_img;i++)
    {
        frecventaB[(liniarizare+i)->B]++;
        frecventaG[(liniarizare+i)->G]++;
        frecventaR[(liniarizare+i)->R]++;
    }
    double chiB=0,chiG=0,chiR=0;
    double chi_mediu=inaltime_img*latime_img;
    chi_mediu=(double)chi_mediu/256;
    for(i=0;i<=255;i++)
    {
        chiB=chiB+(double)((frecventaB[i]-chi_mediu)*(frecventaB[i]-chi_mediu))/chi_mediu;
        chiG=chiG+(double)((frecventaG[i]-chi_mediu)*(frecventaG[i]-chi_mediu))/chi_mediu;
        chiR=chiR+(double)((frecventaR[i]-chi_mediu)*(frecventaR[i]-chi_mediu))/chi_mediu;
    }

    printf("blue:%.2f\ngreen:%.2f\nred:%.2f\n",chiB,chiG,chiR);
    free(frecventaB);
    free(frecventaG);
    free(frecventaR);

}

void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
    FILE *fin, *fout;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char pRGB[3], header[54], aux;

    //printf("nume_fisier_sursa = %s \n",nume_fisier_sursa);

    fin = fopen(nume_fisier_sursa, "rb");
    if(fin == NULL)
    {
        printf("nu am gasit imaginea sursa din care citesc");
        return;
    }

    fout = fopen(nume_fisier_destinatie, "wb+");

    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
    //printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);

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

    //printf("padding = %d \n",padding);

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

double medie_val(unsigned int i,unsigned int j,unsigned int h,unsigned int w,unsigned int hsab,unsigned wsab,pixel*vec)
{
    unsigned int cnt_l,cnt_c;
    double s=0;
    for(cnt_l=i;cnt_l<i+hsab;cnt_l++)
        for(cnt_c=j;cnt_c<j+wsab;cnt_c++)
        {
            s=s+vec[cnt_l*w+cnt_c].R;
        }
    s=s/(wsab*hsab);

    return s;
}


double deviatie_st(unsigned int i,unsigned int j,unsigned int h,unsigned int w,unsigned int hsab,unsigned int wsab,pixel*vec)
{
    double mediu=medie_val(i,j,w,h,wsab,hsab,vec);
    int cnt_l,cnt_c;
    double s=0;
    for(cnt_l=i;cnt_l<i+hsab;cnt_l++)
        for(cnt_c=j;cnt_c<j+wsab;cnt_c++)
        {
            s=s+((vec[cnt_l*w+cnt_c].R-mediu)*(vec[cnt_l*w+cnt_c].R-mediu));
        }
    s=s/((wsab*hsab)-1);
    s=sqrt(s);

    return s;
}


void contur(pixel**I,unsigned int cnt,corelatie*D,unsigned char pR,unsigned char pG,unsigned char pB,unsigned int latime_img,unsigned int inaltime_img,unsigned int latime_sablon,unsigned int inaltime_sablon)
{


    unsigned int k;

    for(k=0;k<latime_sablon;k++)
    {

        (*I+(D[cnt].linie*latime_img+D[cnt].col+k))->R=pR;
        (*I+(D[cnt].linie*latime_img+D[cnt].col+k))->G=pG;
        (*I+(D[cnt].linie*latime_img+D[cnt].col+k))->B=pB;

        (*I+((D[cnt].linie+inaltime_sablon-1)*latime_img+D[cnt].col+k))->R=pR;
        (*I+((D[cnt].linie+inaltime_sablon-1)*latime_img+D[cnt].col+k))->G=pG;
        (*I+((D[cnt].linie+inaltime_sablon-1)*latime_img+D[cnt].col+k))->B=pB;
    }
    for(k=0;k<inaltime_sablon;k++)
    {
        (*I+((D[cnt].linie+k)*latime_img+D[cnt].col))->R=pR;
        (*I+((D[cnt].linie+k)*latime_img+D[cnt].col))->G=pG;
        (*I+((D[cnt].linie+k)*latime_img+D[cnt].col))->B=pB;

        (*I+((D[cnt].linie+k)*latime_img+(D[cnt].col+latime_sablon-1)))->R=pR;
        (*I+((D[cnt].linie+k)*latime_img+(D[cnt].col+latime_sablon-1)))->G=pG;
        (*I+((D[cnt].linie+k)*latime_img+(D[cnt].col+latime_sablon-1)))->B=pB;
    }




}

void ferestre(char*nume_imagine,char*nume_sablon,corelatie**D,int*cnt_ferestre,float ps,unsigned int cifra_sablon,unsigned int *wsab,unsigned int *hsab)
{

    pixel*I,*sablon;
    unsigned char*header;
    citire_imagine(nume_imagine,&I,&header);
    free(header);
    citire_imagine(nume_sablon,&sablon,&header);
    free(header);
    FILE*fimg;
    FILE*fsablon;


    fimg = fopen(nume_imagine, "rb");
    if(fimg == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }


    fsablon = fopen(nume_sablon, "rb");
    if(fsablon == NULL)
    {
        printf("Nu am gasit imaginea sursa(sablonul) din care citesc.");
        return;
    }
    unsigned int dim_img, latime_img, inaltime_img;


    fseek(fimg, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fimg);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fimg, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fimg);
    fread(&inaltime_img, sizeof(unsigned int), 1, fimg);

    unsigned int dim_sablon, latime_sablon, inaltime_sablon;


    fseek(fsablon, 2, SEEK_SET);
    fread(&dim_sablon, sizeof(unsigned int), 1, fsablon);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_sablon);

    fseek(fsablon, 18, SEEK_SET);
    fread(&latime_sablon, sizeof(unsigned int), 1, fsablon);
    fread(&inaltime_sablon, sizeof(unsigned int), 1, fsablon);

    fclose(fimg);
    fclose(fsablon);

    unsigned int i,j,i1,j1;
    double s=0;
    for(i=0;i<inaltime_img&&i+inaltime_sablon<inaltime_img;i++)
        for(j=0;j<latime_img&&j+latime_sablon<latime_img;j++)
        {
            double smediu,imediu,devs,devi;
            imediu=medie_val(i,j,inaltime_img,latime_img,inaltime_sablon,latime_sablon,I);
            smediu=medie_val(0,0,inaltime_sablon,latime_sablon,inaltime_sablon,latime_sablon,sablon);
            devi=deviatie_st(i,j,inaltime_img,latime_img,inaltime_sablon,latime_sablon,I);
            devs=deviatie_st(0,0,inaltime_sablon,latime_sablon,inaltime_sablon,latime_sablon,sablon);

            for(i1=i;i1<i+inaltime_sablon;i1++)
                for(j1=j;j1<j+latime_sablon;j1++)
                {
                    s=s+((I[i1*latime_img+j1].R-imediu)*(sablon[(i1-i)*latime_sablon+(j1-j)].R-smediu));

                }

            s=s/(devi*devs);
            s=s/(inaltime_sablon*latime_sablon);
            if(s>=ps)
            {

                *D=(corelatie*)realloc((*D),(*cnt_ferestre+2)*sizeof(corelatie));
                *cnt_ferestre+=1;
                (*D+(*cnt_ferestre))->linie=i;
                (*D+(*cnt_ferestre))->col=j;
                (*D+(*cnt_ferestre))->corr=s;
                (*D+(*cnt_ferestre))->cifra=cifra_sablon;

            }
            s=0;

        }
    *wsab=latime_sablon;
    *hsab=inaltime_sablon;


}
int compar_corelatie(const void*a,const void*b)
{
     double x=((corelatie*)a)->corr;
     double y=((corelatie*)b)->corr;

     if(x<y)
     return 1;
      else return - 1;
}

void sortare(corelatie**D,int cnt_ferestre)
{
    qsort(*D,cnt_ferestre+1,sizeof(corelatie),compar_corelatie);
}

double arie_intersectie(int i,int j,corelatie*D,unsigned int w,unsigned int h)
{
    double h_arie,w_arie;
    if((abs(D[i].linie-D[j].linie)>=h)||(abs(D[i].col-D[j].col)>=w))
      return 0;
    if(D[i].linie>=D[j].linie)
     h_arie=(D[j].linie+h)-D[i].linie;
       else h_arie=(D[i].linie+h)-D[j].linie;
    if(D[i].col>=D[j].col)
     w_arie=(D[j].col+w)-D[i].col;
       else w_arie=(D[i].col+w)-D[j].col;

   return(h_arie*w_arie);

}

void suprapunere(corelatie**D,int *cnt_ferestre,unsigned int w,unsigned int h)
{
   int i,j,k=0;
   corelatie*D_non_maxime=(corelatie*)malloc(*cnt_ferestre*sizeof(corelatie));
   double arie_int,arie_reun;
   double sup;
   for(i=0;i<=*cnt_ferestre;i++)
    D_non_maxime[i]=(*D)[i];
   for(i=0;i<*cnt_ferestre;i++)
   {
       for(j=i+1;j<=*cnt_ferestre;j++)
       {
           arie_int=arie_intersectie(i,j,*D,w,h);
           arie_reun=2*w*h-arie_int;
           sup=arie_int/arie_reun;

           if(sup>0.2)
           {

               D_non_maxime[j].corr=0;
           }

      }

   }
   k=-1;

   for(i=0;i<=*cnt_ferestre;i++)
   {
       if(D_non_maxime[i].corr!=0)
        {
            k++;
            (*D)[k]=D_non_maxime[i];
        }
   }
   *cnt_ferestre=k;
}

pixel* template_matching(char*nume_imagine,char*nume_grayscale)
{
    unsigned int i;
    int cnt_ferestre=-1;
    pixel*I;
    unsigned char*header;
    unsigned int dim_img, latime_img, inaltime_img;
    grayscale_image(nume_imagine,nume_grayscale);

    FILE*fimg = fopen(nume_imagine, "rb");
    if(fimg == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc.");
        return;
    }

    fseek(fimg, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fimg);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fimg, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fimg);
    fread(&inaltime_img, sizeof(unsigned int), 1, fimg);

    fclose(fimg);

    citire_imagine(nume_imagine,&I,&header);
    free(header);
    corelatie*D=NULL;
    unsigned int wsab,hsab;
    for(i=0;i<10;i++)
    {
        char nume_sablon[40],nume_sablon_gray[40];
        printf("Intrdouceti numele sablonului corespunzator cifrei %u :",i);
        scanf("%s",nume_sablon);
        printf("Introduceti numele imaginii in care va fi salvat sablonul dupa operatia grayscale : ");
        scanf("%s",nume_sablon_gray);
        grayscale_image(nume_sablon,nume_sablon_gray);
        ferestre(nume_grayscale,nume_sablon_gray,&D,&cnt_ferestre,0.5,i,&wsab,&hsab);

    }

    sortare(&D,cnt_ferestre);
    suprapunere(&D,&cnt_ferestre,wsab,hsab);


    for(i=0;i<=cnt_ferestre;i++)
    {

        unsigned int x=D[i].cifra;


        if(x==0) contur(&I,i,D,255,0,0,latime_img,inaltime_img,wsab,hsab);
        if(x==1) contur(&I,i,D,255,255,0,latime_img,inaltime_img,wsab,hsab);
        if(x==2) contur(&I,i,D,0,255,0,latime_img,inaltime_img,wsab,hsab);
        if(x==3) contur(&I,i,D,0,255,255,latime_img,inaltime_img,wsab,hsab);
        if(x==4) contur(&I,i,D,255,0,255,latime_img,inaltime_img,wsab,hsab);
        if(x==5) contur(&I,i,D,0,0,255,latime_img,inaltime_img,wsab,hsab);
        if(x==6) contur(&I,i,D,192,192,192,latime_img,inaltime_img,wsab,hsab);
        if(x==7) contur(&I,i,D,255,140,0,latime_img,inaltime_img,wsab,hsab);
        if(x==8) contur(&I,i,D,128,0,128,latime_img,inaltime_img,wsab,hsab);
        if(x==9) contur(&I,i,D,128,0,0,latime_img,inaltime_img,wsab,hsab);



    }
    return I;

}

int main()
{
    char nume_imagine_initiala[50];
    printf("Introduceti numele imaginii ce urmeaza sa fie criptata : ");
    scanf("%s",nume_imagine_initiala);
    char nume_imagine_criptata[50];
    char cheie_secreta[50];
    char nume_imagine_decriptata[50];
    char nume_test[50];
    char nume_grayscale[50];
    char nume_chenare[50];
    pixel*L,*C,*D;
    pixel*test;
    unsigned char *header;
    citire_imagine(nume_imagine_initiala,&L,&header);
    printf("Introduceti numele fisierului in care se afla cheia secreta : ");
    scanf("%s",cheie_secreta);
    printf("Alegeti un nume pentru imaginea criptata : ");
    scanf("%s",nume_imagine_criptata);
    C=criptare(nume_imagine_initiala,nume_imagine_criptata,cheie_secreta);
    afisare(nume_imagine_initiala,nume_imagine_criptata,C,header);
    printf("Alegeti un nume pentru imaginea decriptata : ");
    scanf("%s",nume_imagine_decriptata);
    D=decriptare(nume_imagine_criptata,nume_imagine_decriptata,cheie_secreta);
    citire_imagine(nume_imagine_criptata,&L,&header);
    afisare(nume_imagine_criptata,nume_imagine_decriptata,D,header);
    printf("Valorile testului chi pentru imaginea initiala : \n");
    chi(nume_imagine_initiala,cheie_secreta);
    printf("\n");
    printf("Valorile testului chi pentru imaginea criptata : \n");
    chi(nume_imagine_criptata,cheie_secreta);
    printf("Introduceti numele imaginii careia i se va aplica functia template matching : ");
    scanf("%s",nume_test);
    printf("Alegeti un nume pentru imaginea grayscale : ");
    scanf("%s",nume_grayscale);
    test=template_matching(nume_test,nume_grayscale);
    citire_imagine(nume_grayscale,&L,&header);
    printf("Alegeti un nume pentru imaginea modificata(dupa operatia template matching si eliminarea non-maximelor):");
    scanf("%s",nume_chenare);
    afisare(nume_grayscale,nume_chenare,test,header);
    free(L);
    free(C);
    free(D);
    free(test);
    return 0;
}

