//
//  PhElem.hpp
//  SDG
//  Elementos fisicos PhElem definido por template
//  para poder especificar o numero de variaveis no momento
//  da sua definicao. 
//
//  Created by Murilo Almeida on 24/10/16.
//  Copyright © 2016 Murilo Almeida. All rights reserved.
//

#ifndef PhElem_hpp
#define PhElem_hpp

#include "spectral.h"

// *************************************************************************
template < int NumVariaveis  >
class PhElem
{
public:
  
  PhElem(); //(const int n = 1);
  ~PhElem(){/* cout << "destruir PhElem\n"; */ };
  //int show_NumLocalVars();
  //void set_NumLocalVars(const int & n);
  void teste();
  void set_ptr_stdel(Stdel * pointer, Stdel * pointer1);
  // void set_ptr_stdel(Mat1<Stdel*> & pointers);
  void set_ptr_stdel(Stdel * pointers[]);
  void set_ptr_stdel(Stdel * pointer);
  void set_ptr_stdel_var(const int ind, Stdel * pointer);
  void set_Vert_map(const int & n_in,int ver_temp[]);
  void set_ptvert(const Vertice * pointer);
  void read_vertices(FILE *, const int & );
  void set_sgn();
  int show_sgn(const int & ia, const int & i){return (sgn[ia][i]);};
  const int show_gbnmap(int ivar, int modo) const { return gbnmap[ivar][modo];};
  int show_numv(){return(numv);};
  int show_nume(){return(nume);};
  int show_numf(){return(numf);};
  void set_type(const int & num);
  int type_val();
  void print_numeracao(FILE * fout,const int & k);
  void print_modes(FILE * fout,const int & k) const;
  void inicia_vetores();
  void inicia_tracos(EDGE *  border);
  void inicia_funcoes_na_borda( EDGE *  border);
  void finaliza_vetores();
  void zera_vetores(const int & k);
  void zera_bs(const int & k);
  void make_vector(const int & ia,double (*f)(double,double,double));
  void make_vector_Elast2D_dyn(const double a[]);
  void make_Kcomp_Elast2D(const double lambda, const double mu);
  void make_Kcomp_Elast2D_dyn(const double lambda,const double mu,
                              const double a[]);
  void make_K_Elast2D(const double lambda, const double mu,double ** K);
  void VectorElast2D(const int & i,double sum[]);
  void VectorElast2D(double vector[]);
  //void AlocarMatrizesK(Matrix * ptKcomp,Matrix * ptKi_inv,Matrix *ptKcKi_inv);
  //void AlocarMatrizesM(Matrix * pt1,Matrix * pt2,Matrix *pt3);
  void P_eval_print_Elast2D(FILE * fout,
                            const double X[],double lambda, double mu);
  void P_eval_print_Elast2D_dyn(FILE * fout,const int & prnflag,
                                const double X[],double lambda, double mu,
                                const double a[]);
  void make_K_Imiscivel2F(const double,const double,const double, const double,
                          const double [],double (*)(double),
                          double (*)(double),double (*)(double));
  double Kcomp(const int & i, const int & j);
  double K_Imisc2F(const int & i, const int & j);
  void make_Vector_Eq_S1(const double, const double,
                         const double, const double,
                         const double [3],double (*)(double));
  //int map(const int & i);
  int map(const int & var,const int & no_local);
  void mapa_inverso(const int & ivar, const double X[]);
  //double MassMatrix_value(int i,int j);
  void P_eval_u(const double X[]);
  void P_eval_u(const double X[],const int & ia);
  void P_eval_u_p1(const double X[],const int & ia);
  void P_eval_print(const double X[],const int & numf,FILE *fout,
                    double (*)(double,double,double));
  void P_print(const int & numf, FILE *file);
  void P_eval_vert(const double X[],double f_vert[]);
  void P_eval_phys(const int & numf,double f0[]);
  void Increment_field(const double X[],const int & numf, const double dt);
  void compute_JV(const int & ia);
  void print_matrices(FILE * fout);
  /*
   void Processar_dados(const int & nel,int& NG,int& NL,std::vector<EDGE>& border,int Ng[],
   int & NF, int Face[], int Fng[], int f_mask[]);
   */
  void Processar_dados(int& NL, std::vector<ARESTA>& aresta,
                       int & NF, std::vector<FACE> & face_vec);//int Face[]);//,int f_mask[]);
  void check_connectivity(FILE * fout,const int & ia);
  void assign_gbnmap(const int & ia,const int & i, const int & val);
  void inicia_gbnmap(int & count);
  void inicia_gbnmap(const int & ivar,int & count);
  void inicia_gbtrbmap(int & count);
  void set_gbnmap(const int & ia,const int gbnmap_temp[],
                  const int sgn_temp[]);
  void set_stgbtrbmap(const int & b, const int & vmapM, const int & vmapP);
  int show_gbtrbmap(){ return(gbtrbmap); };
  int show_stgbtrbmapM(const int & a){ return(stgbtrbmapM[a]); };
  int show_stgbtrbmapP(const int & a){ return(stgbtrbmapP[a]); };
  void check_gradiente(FILE * fout,const int & ia);
  Stdel * show_ptr_stdel(const int & n) const;
  
  //void make_bflag(int bflag[], int facenum);
  void BoundCond(const int & face, int bflag[], double bc[],
                 double (*func)(double,double,double),
                 const int & nf, const int & bndtype);
  void Neumann(const int & aresta,const int & varn,
               const double rho1,const double mu1,
               const double rho2,const double mu2,
               const double gravidade[]);
  double show_u(const int & i);
  double show_u(const int &,const int &);
  double show_b0(const int & i);
  double show_b0(const int &,const int &);
  double show_bs(const int &,const int &);
  void set_mass_density(double);
  void set_porosidade(double);
  void set_fontes(double sw,double sn){qw=sw; qn=sn;};
  
  double show_rho(){return rho;};
  double show_J(){return J;};
  void set_permeabilidade(const double,const double, const double);
  void make_MassMatrices();
  int show_Vert_map(int i){return Vert_map[i];};
  int show_Aresta(int i){return aresta_map[i];};
  int show_Face(int i){return face_map[i];};
  int show_sinal(int i){return sinal[i];};
  int show_border_num(int i){return border_num[i];};
  int show_part_num(){return part_num;};
  double show_Volume(){return Volume;};
  void set_part_num(const int & num = -1);
  double Phi_val(const int & var, const int & ind, const int & pos)
  {return ptr_stdel[var]->show_Phi_val(ind,pos);};
  void projetar_C0(FILE *file,double (*func)(double,double,double),
                   const int & ivar);
  void transformacao_direta(double f[],const int & ivar);
  void teste_transformacao_direta(FILE * fin, FILE * fout, const int & npoints,const double coord[]);
  
  void set_border_bc(/*std::vector<EDGE>&*/EDGE * border,const int & n,const int & t);
  void set_border_num(const int  & aresta, const int & num);
  void Atualizar_valores(FILE * fout = NULL);// NULL eh o valor default
  void Imprimir_valores(FILE * fileout, const int & npoints, const double coord[]);
  void printwGQJ(FILE * fileout);
  void Avancar_u0(const double X[],const double relax = 1.0);
  void Atualizar_u0(const double X[]);
  void Copia_u0_em_(double X[]);
  void Copia_u0_em_(Teuchos::RCP<Epetra_Vector> X);
  void Comparar_u0(const double X[]);
  void Salvar_u0();
  void Restaurar_u0();
  void escrever_restart(FILE * fout);
  void ler_restart(FILE * fin);
  void ler_restart_buffer(FILE * fin,double * buff, int & count);
  void restart_element(double * buff,int & count);
  //void testar_traco(FILE * f_eco);
  void echo_traco(FILE * f_eco = nullptr);
  void teste_gradiente();
  // void teste_tracoI(Fluids fls);
  int  get_trace_border_map(const int & aresta);
  // ***********************************************
  void gbnmap_vertices(std::vector< std::vector<int>> gbn);//int ** gbn);
  void gbnmap_vertices(int * gbn,const int & ivar);
  void gbnmap_faces(const int face_vec[],const int & ivar);
  void gbnmap_aresta(const int & aresta,const int & sinal,const int & ivar,
                     const int & count,int &inc);
  void gbnmap_arestas(const int aresta_vec[],const int & ivar);// nova implementacao
  void gbnmap_interior(const int & ivar,int &count);
  void anexa_gbnmap(const int & ivar, vector <int> &);
  
  void vetor_superficie(const int & num_local,double & area, double normal[3]);
  
  // Funcoes especificas para DG_Problem
  void VolumeIntegrals_UMFPACK(const double Dt,Fluids fls,
                               int & count,
                               int * Ti,int * Tj, double * Tx,
                               double * B,
                               double * = NULL,
                               double * = NULL);
  void VolumeIntegrals(const double Dt,Fluids fls,
                       Teuchos::RCP<Epetra_FECrsMatrix>   A,
                       Teuchos::RCP<Epetra_FEVector> RHS,
                       double * = NULL,
                       double * = NULL);
  void VolumeTracos(const double Dt,Fluids fls,
                    double * = NULL,
                    double * = NULL);
  void VolumeIntegralsT(const double Dt_new,const double Dt_old,int & count,
                        int * Ti,int * Tj, double * Tx,
                        double * B);
  
  void calcula_tracos(Fluids fls);
  // Initial Guess Volume Integrals
  void VolumeIntegrals_IG(Fluids fls, int & count,
                          int * Ti,int * Tj, double * Tx,
                          double * B);
  
  void perm_val(double K[3]);
  
  void Traco_sn(const int & h,double *saida);
  
  void Traco_pw(const int & h,double *saida);
  
  void Traco_Kgrad_pw(const int & h,double ** saida);
  
  void Traco_Kgrad_pc(const int & h,double ** saida);
  
  void Traco_Kgrad_sn(const int & h,double ** saida);
  
  void Traco_phi(const int & lado,const int & ivar,
                 const int & ind,double * saida);
  // ***********************************************************************
  void Traco_phi_1(const int & lado,const int & ivar,const int & pos,
                   double * saida)
  {
    int qmax=ptr_stdel[0]->qborder_val();
    for(int q=0;q<qmax;q++)
      saida[q]=TrPhi[ivar][lado][pos][q];
  };
  // ***********************************************************************
  void Traco_grad_phi(const int & lado,const int & ivar,
                      const int & ind,double ** saida);
  void Traco_Kgrad_phi_n(const int & lado,const int & ivar,
                         const int & ind,double * saida);
  
  // void Row(int * NumNz, int ** MapRow); // Nao esta em uso
  // ***********************************************
  // Versao paralela. Dados em vetor
  // ***********************************************
  // void vetorizar_dados(int elm_num, int &count1, int &count2, int *controle, int *conteudo);
  
  //void construir_PhiArray();
  //const double * eval_Phi(const int i);
  
  
private:
  
  Stdel * ptr_stdel[NumVariaveis];
  int type;
  const Vertice * ptvert;
  int NumLocalVars = NumVariaveis;
  int numv; //!< Number of vertices
  int nume; //!< Number of edges
  int numf; //!< Number of faces
  int numborders; //!< Number of borders in the element
  int Vert_map[8], aresta_map[12],face_map[6], border_num[12], sinal[12];//!< sinal refers to the borders
  int part_num; //!< numero da particao a qual pertence
  int gbnmap[NumVariaveis][MAXNN]; //!< mapping from local to global nodes (or modes)
  int sgn[NumVariaveis][MAXNN];  //!< sign of the local mode (or mode)
  double J; //!< Jacobian from the physical to the standard element
  double * JV; //!<Ponteiro para matriz do Jacobiano nos pontos gaussianos;
  double * b0[NumVariaveis];     //!<[MAXMODES*MAXNFIELDS];  vector of body forces
  double * bs[NumVariaveis];     //!< vector of force boundary conditions
  double * u0[NumVariaveis]; //!<vector with the coefficients of each mode
  double * usave[NumVariaveis];
  //double * ua[NumVariaveis];
  
  int gbtrbmap ;/*!< \brief Beginning index on the global array where
                 * the local border Gauss quadrature points (traces)
                 * are mapped; this map facilitates
                 * the use of mpi for paralell calculations over
                 * element borders */
  int * stgbtrbmapM; /*!< \brief Array containing the start point in the global
                      *  trace array of the internal trace of the border */
  int * stgbtrbmapP; /*!< \brief Array containing the start point in the global
                      *  trace array of the external trace of the border */
  
  // Variaveis especificas para o caso DG_Problem
  double * sna, * pwa; //!< ponteiros para os valores nos pontos de Gauss da saturacao e pressao atuais
  double rho; //!< densidade de massa ou porosidade
  double porosidade;
  double perm[3]; //!< permeabilidades nas 3 direcoes
  double qn,qw;
  double Volume; //!< Volume do elemento ( = area do elemento para elementos bidimensionais e comprimento para elem. 1D
  int vetores_iniciados; // = 1; indica que os vetores locais foram iniciados e necessitam ser finalizados
  double     * Jb; //!< Jacobiano nos pontos de Gauss sobre as bordas; usado na integração sobre bordas
  double    ** Mass_sn;
  double    ** Trsn;
  double    ** Trpw;
  double   *** LaplacianoPhi;
  double   *** TrKgrad_sn;
  double   *** TrKgrad_pw;
  double   *** TrKgrad_pc;
  double  **** TrPhi;
  double  **** GradPhi;
  double  **** TrKgradPhi_n;
  double ***** TrGradPhi;
  double     * PhiArray;
  
  //double ** b;
  
};
/*! \class PhElem
 * Physical Elements
 */


// ****************************************************************************
// Class PhElem ---  Para escoamentos de 2 fluidos imisciveis DG
// ****************************************************************************

template < int NumVariaveis >
PhElem<NumVariaveis>::PhElem()
{
  //set_ptr_stdel(NULL,NULL);
  rho=1.0;
  vetores_iniciados=0;// 0 (=FALSE)
  // qn=0.0;
  // qw=0.0;
  // default value
  //printf("Inicializa PhElem: NumLocalVars = %d rho = %g\n",n,rho);
};
// ****************************************************************************

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::escrever_restart(FILE * fout)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) fprintf(fout," %e",u0[ivar][i]);
  }
  fprintf(fout,"\n");
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::ler_restart(FILE * filein)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) fscanf(filein,"%lf",&u0[ivar][i]) ;
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::ler_restart_buffer(FILE * filein,double * buff, int & contador)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) fscanf(filein,"%lf",&buff[contador++]) ;
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::restart_element(double * buff,int & conta)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) u0[ivar][i] = buff[conta++] ;
  }
};


// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::mapa_inverso(const int & ivar, const double X[])
{
  double temp;
  int nn=ptr_stdel[ivar]->nn_val();
  for(int i=0;i<nn;i++){
    temp = X[gbnmap[ivar][i]];
    u0[ivar][i]= temp;
    if(std::isnan(temp))cout << "Mapa inverso Float was Not a Number: ivar " << ivar << " modo "<< i << endl;
  }
  
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Copia_u0_em_(double X[])
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) X[gbnmap[ivar][i]] = u0[ivar][i];
  }
};
template<int NumVariaveis>
void PhElem<NumVariaveis>::Copia_u0_em_(Teuchos::RCP<Epetra_Vector> X)
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) (*X)[gbnmap[ivar][i]] = u0[ivar][i];
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::anexa_gbnmap(const int & ivar, vector <int> & list)
{
  int nn=ptr_stdel[ivar]->nn_val();
  for(int i=0;i<nn;i++) list.push_back( gbnmap[ivar][i] );
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Atualizar_valores(FILE * fileout)
{
  
  //cout << "sat " << endl;
  ptr_stdel[sat]->evalGQ(sna,u0[sat]);
  //cout << "pres " << endl;
  ptr_stdel[pres]->evalGQ(pwa,u0[pres]);
  
  if(fileout != NULL){
    ptr_stdel[sat]->printGQtofile(fileout,sna,pwa,ptvert,Vert_map);
    fflush(fileout);
  }
};

// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::printwGQJ(FILE * fileout)
{
  if(fileout != NULL){
    ptr_stdel[0]->printwGQtofile(fileout,ptvert,Vert_map,JV);
    fflush(fileout);
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Salvar_u0()
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) usave[ivar][i] = u0[ivar][i];
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Restaurar_u0()
{
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    for(int i=0;i<nn;i++) u0[ivar][i] = usave[ivar][i];
  }
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Avancar_u0(const double X[],const double relax)
{
  // default relax = 1.0
  
  // atualizar os coeficientes: u0 -= X
  
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    // *****************************************************************
    for(int i=0;i<nn;i++){                                          // *
      u0[ivar][i] -= relax * X[gbnmap[ivar][i]];                    // *
    }                                                               // *
    // *****************************************************************
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Comparar_u0(const double X[])
{
  // default relax = 1.0
  
  
  
  // comparar os coeficientes: u0 e X
  
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    // *****************************************************************
    for(int i=0;i<nn;i++){                                          // *
      if(u0[ivar][i] != X[gbnmap[ivar][i]])
        cout << "Falhou em VarGlobal " << gbnmap[ivar][i] << "\n" ; // *
    }                                                               // *
    // *****************************************************************
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Atualizar_u0(const double X[])
{
  
  // atualizar os coeficientes: u0 = X
  
  for(int ivar=0;ivar<NumVariaveis;ivar++){
    int nn=ptr_stdel[ivar]->nn_val();
    // *****************************************************************
    for(int i=0;i<nn;i++){                                          // *
      u0[ivar][i] = X[gbnmap[ivar][i]];                             // *
    }                                                               // *
    // *****************************************************************
  }
};
/*
 // ****************************************************************************
 template<int NumVariaveis>
 void  PhElem::testar_traco(FILE * f_eco)
 {
 const int qmax=ptr_stdel[0]->qborder_val();
 const int ndim=ptr_stdel[0]->ndim_val();
 int NGQP=ptr_stdel[0]->NGQP_val();
 double * gphi_r[ndim];
 double * gphi_l[ndim];
	double A[2][qmax];
	double B[2][qmax];
 for(int i=0; i<2; i++){
 gphi_r[i]=A[i];
 gphi_l[i]=B[i];
 }
 
 fprintf(f_eco,"Phi para o elemento\n");
 for(int i=0;i<2;i++){
 fprintf(f_eco,"\n variavel %d\n",i);
 int nn=ptr_stdel[i]->nn_val();
 int nborder=ptr_stdel[i]->nborder_val();
 for(int j=0;j<nn;j++){
 fprintf(f_eco,"         modo %d\n",j);
 fprintf(f_eco,"Grad_Phi[%d][%d]\n",i,j);
 for(int m=0;m<NGQP;m++){
 fprintf(f_eco," m= %d ",m);
 for (int ndir=0; ndir< ndim; ndir++){
 fprintf(f_eco,"%11.4e ",GradPhi[i][j][ndir][m]);
 }
 fprintf(f_eco,"\n");
 }
 printf("AQUI1 *********************************************************\n");
 for(int h=0;h<nborder;h++){
 Traco_grad_phi(h,i,j,gphi_r);
 fprintf(f_eco,"  Traco na aresta %d \n",h);
 for(int m=0;m<qmax;m++){
 for (int ndir=0; ndir< ndim; ndir++){
 fprintf(f_eco,"%11.4e ",gphi_r[ndir][m]);
 }
 fprintf(f_eco,"\n");
 }
 }
 printf("AQUI2 *******************************************************\n");
 }
 }
 // **********************************************************
 // Fim do echo dos resultados - testar_traco
 // **********************************************************
 };
 */

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_permeabilidade(const double x,const double y, const double z)
{
  perm[0]=x;
  perm[1]=y;
  perm[2]=z;
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_mass_density(double value)
{
  rho=value;
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_porosidade(double value)
{
  porosidade=value;
};
// ****************************************************************************
//template<int NumVariaveis>
//void PhElem<NumVariaveis>::set_NumLocalVars(const int & n)
//{
// NumLocalVars=n;
//};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel(Stdel * point0 )
{
  ptr_stdel[0]=point0;
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel(Stdel * point0, Stdel * point1 )
{
  if (NumVariaveis > 2) {
    ptr_stdel[0]=point0;
    ptr_stdel[1]=point1;
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel_var(const int var, Stdel * point0 )
{
  ptr_stdel[var]=point0;
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptr_stdel(Stdel *  ponteiros[])
{
  for(int i=0;i<NumVariaveis;++i)
    ptr_stdel[i]=ponteiros[i];
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_ptvert(const Vertice * pointer)
{
  ptvert=pointer;
};
// ****************************************************************************
// Especifica os valores das condicoes de contorno
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_border_bc( EDGE * border,const int & aresta,const int & t)
{
  int ndim = ptr_stdel[0]->ndim_val();
  int en = border_num[aresta];
  // printf("entrou PhElem<NumVariaveis>::set_border_bc: border_num=%d\n",en);
  border[en].tipo=t;
  
  double x,y,aux;
  
  switch (ndim) {
      
    case 1:
      
      //printf("Entrando set_border_bc\n");
      if(t == -1 || t == 1){
        int na=border[en].Na;
        x=ptvert[na].x;
        y=0.0;
        border[en].pdir = new double [1];
        border[en].pdir[0]=funcao_pdir(x,y,t);
        if(t==-1){
          border[en].sdir = new double [1];
          border[en].sdir[0]=funcao_sdir(x,y);
        }
      }
      //printf("Saindo de set_border_bc edge number = %d\n",en);
      break;
      
    case 2:
      
      if(t == -1 || t == 1){
        // Alocar memoria para condicoes de contorno de Dirichlet
        int qmax=ptr_stdel[0]->qborder_val();
        double xq[qmax],w[qmax];
        double Dtemp[MAXQ][MAXQ];
        //Mat2<double> Dtemp(qmax,qmax);
        Gauss_Jacobi_parameters(qmax,0.0,0.0,xq,w,Dtemp);
        // *******************************************************
        int na=border[en].Na;
        int nb=border[en].Nb;
        double xa=ptvert[na].x;
        double ya=ptvert[na].y;
        double xb=ptvert[nb].x;
        double yb=ptvert[nb].y;
        double xsum=(xb+xa)*0.5;
        double xdif=(xb-xa)*0.5;
        double ysum=(yb+ya)*0.5;
        double ydif=(yb-ya)*0.5;
        
        border[en].pdir = new double [qmax];
        for(int i=0;i<qmax;i++){
          aux=xq[i];
          x=xsum+xdif*aux;
          y=ysum+ydif*aux;
          border[en].pdir[i]=funcao_pdir(x,y,t);
        }
        if(t==-1){
          border[en].sdir = new double[qmax];
          for(int i=0;i<qmax;i++){
            aux=xq[i];
            x=xsum+xdif*aux;
            y=ysum+ydif*aux;
            border[en].sdir[i]=funcao_sdir(x,y);
          }
        }
      }
      break;
  } // switch (ndim)
  //printf("saindo de PhElem<NumVariaveis>::set_border_bc\n");
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_border_num(const int & aresta, const int & num)
{
  border_num[aresta]=num;
};
// ****************************************************************************
template<int NumVariaveis>
int PhElem<NumVariaveis>::get_trace_border_map(const int & aresta)
{
  int qmax=ptr_stdel[0]->qborder_val();
  int temp = gbtrbmap+qmax*aresta;
  return temp;
}
;
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_type(const int & num)
{
  type = num;
};
// ****************************************************************************
template<int NumVariaveis>
int PhElem<NumVariaveis>::type_val()
{
  return type;
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::print_numeracao(FILE * fout,const int & k)
{
  int nb=ptr_stdel[k]->nb_val(); // number of  boundary modes
  fprintf(fout,"Nos nas arestas da variavel %d\n",k);
  for(int i=0; i<nb;i++){
    fprintf(fout,"print_numeracao nl=%4d ng=%4d\n",i,gbnmap[k][i]);
  }
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::print_modes(FILE * fout,const int & ia) const
{
  int nv=ptr_stdel[ia]->nv_val(); // number of vertices
  fprintf(fout,"modes da variavel %d = %d",ia,nv);
  for(int i=0; i<nv;i++){
    fprintf(fout," %d",gbnmap[ia][i]);
  }
  fprintf(fout,"\n");
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::read_vertices(FILE * finput, const int & nv) // nao usada
{
  int label;
  for(int i=0;i<nv;i++)
    fscanf(finput,"%d",&Vert_map[i]);
  fscanf(finput,"%d",&label);
};

// ****************************************************************************
template<int NumVariaveis>
int PhElem<NumVariaveis>::map(const int & var,const int & no_local)
{
  return (gbnmap[var][no_local]);
};

// ****************************************************************************

//template<int NumVariaveis>
//int PhElem<NumVariaveis>::show_NumLocalVars()
//{
//  return(NumLocalVars);
//};

// ***************************************************************************
template<int NumVariaveis>
Stdel * PhElem<NumVariaveis>::show_ptr_stdel(const int & n) const
{
  return ptr_stdel[n];
};

// ****************************************************************************
template<int NumVariaveis>
double PhElem<NumVariaveis>::show_u(const int & var,const int & no)
{
  return u0[var][no];
};

// ****************************************************************************
template<int NumVariaveis>
double PhElem<NumVariaveis>::show_b0(const int & var, const int & no)
{
  return (b0[var][no]);// * sgn[var][no]);
};

// ****************************************************************************
template<int NumVariaveis>
double PhElem<NumVariaveis>::show_bs(const int & var, const int & no)
{
  return (bs[var][no]);// * sgn[var][no]);
};


// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::P_eval_phys(const int & nvar,double f0[])
{
  int nn=ptr_stdel[nvar]->nn_val();
  int i;
  double utemp[nn];
  for(i=0;i<nn;i++){
    utemp[i]=u0[nvar][i];
  }
  ptr_stdel[nvar]->evalGQ(f0,utemp);
};


// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::P_print(const int & ivar, FILE *file)
{
  int nn=ptr_stdel[ivar]->nn_val();
  int i;
  double utemp[nn];
  for(i=0;i<nn;i++){
    utemp[i]=u0[ivar][i];
  }
  ptr_stdel[ivar]->printtofile(file,utemp,ptvert,Vert_map);
};


// ****************************************************************************
// Processar dados dos elementos apontando para                               *
// a rotina de cada tipo de elemento                                          *
// ****************************************************************************
/* //se tornou obsoleta
 template<int NumVariaveis>
 void  PhElem::Processar_dados(const int & nel,int& NG,int& NL,
 std::vector<EDGE>& border,
 int Ng[], //inicio da numeracao da borda
 int & NF, int Face[], int Fng[], //inicio da numeracao da face
 int f_mask[])
 {
 int gbnmap_temp[MAXNN]; // mapping from local to global modes (or modes)
 int sgn_temp[MAXNN];
 int sinal_temp[12];
 cout << "PhElem::Processar_dados (velha)\n";
 ptr_stdel[0]->Processar_geometria(nel,ptvert,numv,Vert_map,
 gbnmap_temp,sgn_temp, // candidato a sair
 sinal_temp,
 NG,NL,border,
 Ng, // candidato a sair
 NF,Face,
 Fng,  // candidato a sair aguardar ver se necessario em Tetrahedral
 f_mask);
 
 //  for(int i=0;i<ptr_stdel[0]->nn_val(); i++){
 //     gbnmap[0][i]=gbnmap_temp[i]; // gbnmap do modo
 //     sgn[0][i]=sgn_temp[i]; // sinal do modo
 //   }
 
 for(int i=0;i<ptr_stdel[0]->ne_val();i++) sinal[i]=sinal_temp[i]; // sinal da aresta
 //printf("PhElem::Ler PONTO 1c\n");
 };
 */

// **************************************************************************
// Candidato a substituir a funcao acima
//
template<int NumVariaveis>
void PhElem<NumVariaveis>::Processar_dados(int& NL,
                                                 std::vector<ARESTA>&aresta_vec,
                                                 int & NF,
                                                 std::vector<FACE> & face_vec)
{
  int na,nb;
  //cout << "PhElem<NumVariaveis>::Processar_dados (nova)\n";
  // recuperar os numeros de vertices, arestas e faces do elemento padrao
  int ndim = ptr_stdel[0]->ndim_val();
  numv=ptr_stdel[0]->nv_val();
  nume=ptr_stdel[0]->ne_val();
  numf=ptr_stdel[0]->nf_val();
  
  // Processar as arestas
  if(ndim!=1) {
    for(int i = 0; i<nume;++i){
      int n0=Vert_map[ptr_stdel[0]->aresta_lvert(i,0)];
      int n1=Vert_map[ptr_stdel[0]->aresta_lvert(i,1)];
      if(n0<n1){ // vertices em ordem crescente
        sinal[i]=1;
        na=n0;
        nb=n1;
      }
      else { // vertices em ordem decrescente
        sinal[i]=-1;
        na=n1;
        nb=n0;
      }
      aresta_map[i] = aresta_gbnum(na,nb,NL,aresta_vec);
      // cout << "aresta "<< aresta_map[i] << " ("<< na<< ","<< nb<< ")"<< endl;
    }
  }
  
  // Processar as faces
  if(ndim==3) {
    for(int i=0;i<numf;++i){
      int nvf = ptr_stdel[0]->show_nvf(i);
      int var[nvf];
      for(int j=0;j<nvf;++j){
        var[j]= Vert_map[ptr_stdel[0]->face_lvert(i,j)];
      }
      int n = face_gbnum(nvf,var,NF,face_vec);
      face_map[i] = n;
      if(face_vec[n]._tipo == 0){
        face_vec[n]._tipo=ptr_stdel[0]->show_face_tipo(i);
      }
    }
  }
  //printf("PhElem<NumVariaveis>::Ler PONTO 1c\n");
};

// ***************************************************************************
// Salva o numero de vertices e o o mapa de Vertices Vert_map
// ***************************************************************************
// 19/07/2014
/*
 template<int NumVariaveis>
 void PhElem::set_Vert_map(const int & n_in,int ver_temp[])
 {
 numv=n_in;
 
 if(type==4){ // tetrahedral
 int seq_orig[4] = {0,1,2,3};
 ordenar(4,ver_temp,seq_orig);
 face_mask = new int[4];
 int facemap[4][3] = // sequencia das faces do tetraedro
 { 0,1,2,
 0,1,3,
 1,2,3,
 0,2,3
 };
 int temp[3],temp1[3] = {0,1,2};
 // face nova i em termos de seq_orig
 for(int i=0; i<4; ++i) { // nova face
 for(int j=0; j<3; ++j) temp[j] =seq_orig[facemap[i][j]]; //nova face = i
 ordenar(3,temp,temp1);
 // comparar se temp se iguala a face m
 int flag=0;
 for(int m=0; m<4 && flag == 0; ++m){ //face antiga = m
 flag=1;
 for(int k=0; k<3 && flag ==1; ++k){
 if(temp[k] != facemap[m][k]) flag = 0; // faces diferem
 }
 if(flag == 1) face_mask[m]=i; //faces coincidem
 }
 }
 }
 }
 */
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_Vert_map(const int & n_in,int ver_temp[])
{
  numv=n_in;
  for(int i=0;i<numv;++i)Vert_map[i]=ver_temp[i];
  if(type==5) { // hexaedro requer ordenacao dos nos lidos de gmsh;
    Vert_map[2]=ver_temp[3];
    Vert_map[3]=ver_temp[2];
    Vert_map[6]=ver_temp[7];
    Vert_map[7]=ver_temp[6];
  }
}
// ***************************************************************************
// Salva os mapas (local para global) e os sinais
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_gbnmap(const int & ia,const int gbnmap_temp[],
                                            const int sgn_temp[])
{
  for(int i=0;i<ptr_stdel[ia]->nn_val(); i++){
    gbnmap[ia][i]=gbnmap_temp[i];
    sgn[ia][i]=sgn_temp[i];
  }
};

// *********************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::check_connectivity(FILE * fout,const int & ia)
{
  int i,p,q,r;
  for(i=0; i<ptr_stdel[ia]->nb_val(); i++){
    ptr_stdel[ia]->show_ind(i,p,q,r);
    fprintf(fout,"i = %3d gbnmap = %3d p = %3d q = %3d\n", i, gbnmap[ia][i],p,q);
  }
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::check_gradiente(FILE * fileout,const int & ia)
{
  int NGQP=ptr_stdel[ia]->NGQP_val();
  int ndim=ptr_stdel[ia]->ndim_val();
  double *ptr_grad[ndim];
  double grad[ndim][NGQP];
  double sn[NGQP];
  for(int i=0;i<ndim;i++){
    ptr_grad[i]=grad[i];
  }
  ptr_stdel[ia]->evalGQ(sn,u0[ia]);
  ptr_stdel[ia]->Gradiente(fileout,ptr_grad,sn,ptvert,Vert_map);
  
  //ptr_stdel[ia]->Gradiente(fileout,ptr_grad,funcao,ptvert,Vert_map);
  // for(int m=0;m<NGQP;m++)
  //fprintf(fileout,"check_gradiente %11.4e\n",grad[0][m]);
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::assign_gbnmap(const int & ia,const int & i, const int & val)
{
  gbnmap[ia][i]=val;
  //  sgn[ia][i]=1;
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_gbnmap(int & count)
{
  for(int k=0;k<NumVariaveis;++k){
    inicia_gbnmap(k,count);
  }
};
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_gbnmap(const int & ivar,int & count)
{ // Nao impoe continuidade entre os elementos vizinhos
  int NN=ptr_stdel[ivar]->nn_val();
  for(int i=0;i<NN;++i){
    gbnmap[ivar][i]=count++;
    sgn[ivar][i]=1;
  }
};
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_gbtrbmap(int & count)
{
  // qmax= number of quadrature points on each border
  const int N = ptr_stdel[0]->nborder_val();
  const int qmax= ptr_stdel[0]->qborder_val();
  //a ser substituido por numborders(number of borders)
  stgbtrbmapM = new int [N];
  stgbtrbmapP = new int [N];
  gbtrbmap = count;
  count += (N*qmax);
  // printf(" gbtrbmap %d\n", gbtrbmap);
};
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_stgbtrbmap(const int & b,const int & vmapM, const int & vmapP)
{
  stgbtrbmapM[b] = vmapM;
  stgbtrbmapP[b] = vmapP;
};

// ******************************************************************************/
// Pressupoem que o JV (vetor dos Jacobianos nos pontos de Gauss ja foi calculado
// ******************************************************************************/
template<int NumVariaveis>
void PhElem<NumVariaveis>::projetar_C0(FILE *file,
                                             double (*func)(double,double,double),
                                             const int & ivar)
{
  //printf("\nComeco de PhElem::projetar_C0\n");
  int i,j,ii,jj,k;
  int nborder=ptr_stdel[ivar]->nborder_val();
  int nb=ptr_stdel[ivar]->nb_val();
  int nn=ptr_stdel[ivar]->nn_val();
  int ni=nn-nb;
  //cout << "PhElem::projetar_C0 ni = "<< ni << endl;
  double aux=0.0;
  int nmap[nn],bflag[nn],sgn[nn];
  double Xl[nn];// coeficientes dos modos locais; vetor a ser calculado
  //cout << "\n Projetar C0\n";
  // cout << "nb "<< nb << " nn "<< nn <<" nborder "<< nborder<< "\n";
  for(i=0;i<nn;i++){
    nmap[i]=i;
    sgn[i]=1;
    bflag[i]=1;
    Xl[i]=0.0;
  }
  
  // for(int i=0;i<ptr_stdel[0]->nv_val();++i)
  //   cout << "Vert_map["<< i<< "] = "<<Vert_map[i]<<endl;
  for(j=0;j<nborder;j++){
    //j=3;
    //cout << "Chamar Dirichlet (ptr_stdel) para border = "<< j << endl;
    
    ptr_stdel[ivar]->Dirichlet(j,ptvert,Vert_map,nmap,sgn,bflag,Xl,func);
    
    // cout << "saindo de Dirichlet para a face "<< j << endl;
    //for(int k=0;k<nn;++k) cout << "Xl["<< k << "%d] = " <<Xl[k] << endl;
  }
  //cout << "Terminou loop sobre as bordas (ptr_stdel)\nRetornou a PhElem::projetar\n";
  
  //ni=0;
  if(ni>0){
    //cout << "calcular o produto interno de f por phi dos "<< ni << " modos internos"<< endl;
    
#ifdef _NEWMAT
    // newmat
    NEWMAT::Matrix Mi(ni,ni);
    NEWMAT::ColumnVector B(ni);
#endif
    
    //fprintf(file,"\n\n%5d\n",ni*ni);
    //printf("Matriz Mi\n");
    for(i=nb;i<nn;++i){
      ii=i-nb;
      Mi.element(ii,ii)=ptr_stdel[ivar]->mass(i,i,JV);
      // printf("%3d %3d %g\n",ii,ii,Mi.element(ii,ii));
      for(j=i+1;j<nn;j++){
        jj=j-nb;
        aux=ptr_stdel[ivar]->mass(i,j,JV);
        Mi.element(ii,jj)=aux;
        Mi.element(jj,ii)=aux;
        //printf("%3d %3d %g\n%3d %3d %g\n",ii,jj,aux,jj,ii,aux);
      }
    }
    
    int NGQP = ptr_stdel[ivar]->NGQP_val();
    double phi[NGQP];
    double f[NGQP];
    // inicializar f;
    ptr_stdel[ivar]->computeFuncGQ(f,ptvert,Vert_map,func);
    for(j=0;j<nb;j++){
      ptr_stdel[ivar]->eval_Phi(j,phi);
      double aux = Xl[j];
      for(k=0;k<NGQP;k++) {f[k] -= (aux*phi[k]);}
    }
    double b[nn];
    ptr_stdel[ivar]->inner_product_vector(b,f,JV);
    //fprintf(file,"%5d\n",ni);
    for(i=0;i<ni;i++){
      B.element(i)=b[i+nb];
      //fprintf(file,"%g ",B[i]);
    }
    //fprintf(file,"\n");
    
#ifdef _NEWMAT
    NEWMAT::ColumnVector Y = Mi.i() * B;  B=Y ; // newmat
#endif
    
    for(i=nb;i<nn;++i)u0[ivar][i]=B.element(i-nb);
  }
  // fim de if(ni>0)
  
  for(i=0;i<nb;++i) u0[ivar][i]=Xl[i]; // no contorno
  
  if(file != NULL)
    ptr_stdel[ivar]->printtofile(file,u0[ivar],func,ptvert,Vert_map);
  
  //cout << "Terminou  PhElem::projetar_C0\n\n";
  
};

// ***********************************************************************
// calcula os coeficientes a partir dos valores nos pontos de Gauss      *
// ***********************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::transformacao_direta(double f[],const int & ivar)
{
  int i,j;
  int nn=ptr_stdel[ivar]->nn_val();
  double aux;
  //cout << "\n Projetar C0\n";
  
#ifdef _NEWMAT
  // newmat
  NEWMAT::Matrix Mi(nn,nn);
  NEWMAT::ColumnVector B(nn);
#endif
  
  //fprintf(file,"\n\n%5d\n",ni*ni);
  for(i=0;i<nn;i++){
    Mi.element(i,i)=ptr_stdel[ivar]->mass(i,i,JV);
    //fprintf(file,"%3d %3d %g\n",ii,ii,Mi[ii][ii]);
    for(j=i+1;j<nn;j++){
      aux=ptr_stdel[ivar]->mass(i,j,JV);
      Mi.element(i,j)=aux;
      Mi.element(j,i)=aux;
      //fprintf(file,"%3d %3d %g\n%3d %3d %g\n",ii,jj,aux,jj,ii,aux);
    }
  }
  
  double b[nn];
  ptr_stdel[ivar]->inner_product_vector(b,f,JV);
  for(i=0;i<nn;i++){
    B.element(i)=b[i];
  }
#ifdef _NEWMAT
  NEWMAT::ColumnVector Y = Mi.i() * B;  B=Y ; // newmat
#endif
  
  for(i=0;i<nn;i++)u0[ivar][i]=B.element(i);
};

// ****************************************************************************
// template<int NumVariaveis>
// void PhElem<NumVariaveis>::BoundCond(const int & face,int bflag[],double X[],
// 		     double (*func)(double,double,double),
// 		     const int & nf, const int & bndtype)
// {
//   // ***********************************************************************
//   // bndtype = 0 Dirichlet; marca bflag = 0 (variavel conhecida)
//   // bndtype = 1 Especifica forca; marca bflag = 1 (variavel desconhecida)
//   // ***********************************************************************
//   if(bndtype==0)
//     ptr_stdel[nf]->Dirichlet(face,ptvert,gbnmap[nf],sgn[nf],bflag,X,func,nf,NumLocalVars);
//   if(bndtype==1){
//     printf("BNDTYPE = 1\n");
//     ptr_stdel[nf]->BoundForce(face,ptvert,gbnmap[nf],bs[nf],func,nf,NumLocalVars);
//     }
// }


//  // ************************************************************************
//  // ************************************************************************
//  template<int NumVariaveis>
//  void PhElem<NumVariaveis>::AlocarMatrizesK(Matrix * pt1,Matrix * pt2,Matrix *pt3)
//  {
//    ptKcomp=pt1;
//    ptKi_inv=pt2;
//    ptKcKi_inv=pt3;
//  };

//    // **********************************************************************
//    template<int NumVariaveis>
//    void PhElem<NumVariaveis>::AlocarMatrizesM(Matrix * pt1,Matrix * pt2,Matrix *pt3)
//    {
//      ptMb_comp=pt1;
//      ptMi_inv=pt2;
//      ptMcMi_inv=pt3;
//    };
//    // **********************************************************************
//    template<int NumVariaveis>
//    void PhElem<NumVariaveis>::make_MassMatrices()
//    {
//      int i,j,ii,jj,im,jm;
//      int nb=ptr_stdel->nb_val();
//      int Nb=nb*NumLocalVars;
//      int nn=ptr_stdel->nn_val();
//      int Nn=nn*NumLocalVars;
//      int Ni=Nn-Nb;
//      int ni=nn-nb;
//
//      // Alocacao dinamica das matrizes
//      ptMb_comp   = new Matrix(0,Nb-1,0,Nb-1);
//      ptMi_inv  = new Matrix(0,Ni-1,0,Ni-1);
//      ptMcMi_inv= new Matrix(0,Nb-1,0,Ni-1);
//
//      for(int i=0;i<Nb;i++){
//        for(int j=0;j<Nb;j++)
//          (*ptMb_comp)(i,j)=0.0;
//        for(int j=0;j<Ni;j++)
//          (*ptMcMi_inv)(i,j)=0.0;
//      }
//      for(int i=0;i<Ni;i++)
//        for(int j=0;j<Ni;j++)
//          (*ptMi_inv)(i,j)=0.0;
//
//      NEWMAT::Matrix mb(0,nb-1,0,nb-1);
//      NEWMAT::Matrix mb_c(0,nb-1,0,nb-1);
//      NEWMAT::Matrix mi(0,ni-1,0,ni-1);
//      NEWMAT::Matrix mc(0,nb-1,0,ni-1);
//      NEWMAT::Matrix mi_inv(0,ni-1,0,ni-1);
//      NEWMAT::Matrix mcmi_inv(0,nb-1,0,ni-1);
//
//      // Matriz Mb
//      for(i=0;i<nb;i++){
//        for(j=0;j<nb;j++) mb.element(i,j)=ptr_stdel->mass(i,j,JV);
//        // Matrix Mc
//        for(j=nb;j<nn;j++){
//          jj=j-nb;
//          mc.element(i,jj)=ptr_stdel->mass(i,j,JV);
//        }
//      }
//      // Matrix Mi
//      for(i=nb;i<nn;i++){
//        ii=i-nb;
//        mi.element(ii,ii)=ptr_stdel->mass(i,i,JV);
//        for(j=i+1;j<nn;j++){
//          jj=j-nb;
//          double aux=ptr_stdel->mass(i,j,JV);
//          mi.element(ii,jj)=aux;
//          mi.element(jj,ii)=aux;
//        }
//      }
//
//      mi_inv=mi.Inverse();
//      mcmi_inv = mc * mi_inv;
//      mb_c = mb - (mcmi_inv * mc.Transpose());
//
//      for(int ivar=0;ivar<NumLocalVars;ivar++){
//        for(i=0;i<nb;i++){
//          im=i*NumLocalVars+ivar;
//          for(j=0;j<nb;j++)
//    	(*ptMb_comp)[im][j*NumLocalVars+ivar]=mb_c.element(i,j);
//          for(j=nb;j<nn;j++){
//    	jj=j-nb;
//    	jm=jj*NumLocalVars+ivar;
//    	(*ptMcMi_inv)[im][jm]=mcmi_inv.element(i,jj);
//          }
//        }
//        //printf("Making local matrices 3\n");
//        for(i=nb;i<nn;i++){
//          ii=i-nb;
//          im=ii*NumLocalVars+ivar;
//          for(j=nb;j<nn;j++){
//    	jj=j-nb;
//    	jm=jj*NumLocalVars+ivar;
//    	(*ptMi_inv)[im][jm]=mi_inv.element(ii,jj);
//          }
//        }
//      }
//    }


// ****************************************************************************
// PhElem<NumVariaveis>::Vector(int i) creates the vector for the boundary modes performing
// the Schur condensation Vector = fb-McMi_inv*fi and applying the sgn[i]
// ****************************************************************************
// template<int NumVariaveis>
// double PhElem<NumVariaveis>::Vector(const int & ia,const int & i)
// {
//     int nb=ptr_stdel[ia]->nb_val();
//     int ni=(ptr_stdel[ia]->nn_val()-nb);
//     double sum0=b0[ia][i];
//     for(int j=0; j<ni; j++)
//       sum0-=(*ptMcMi_inv)(i,j)*b0[ia][j+nb];
//     return (sum0*sgn[ia][i]);
// };

//  // ***********************************************************************
// template<int NumVariaveis>
//  void PhElem<NumVariaveis>::Vector2(const int & i,double sum[])
//  {
//    // *********************************************************
//    // i = numero local do no; cada no tem NumLocalVars variaveis *
//    // *********************************************************
//    //printf("entrou em vector2\n");
//    int nb=ptr_stdel->nb_val();// numero de nos de contorno
//    int ni=ptr_stdel->nn_val()-nb;// numero de nos internos
//    int Ni=ni*NumLocalVars;
//    int Nb=nb*NumLocalVars;
//    int ia;
//    for(ia=0;ia<NumLocalVars;ia++){
//      int ii = NumLocalVars*i + ia;
//      sum[ia]=b0[ii];
//      for(int j=0; j<Ni; j++){
//        sum[ia]-=(*ptMcMi_inv)(ii,j)*b0[j+Nb];
//      }
//    }
//    for(ia=0;ia<NumLocalVars;ia++)
//      sum[ia]*=sgn[i];
//  };

//  // ************************************************************************
//  double PhElem<NumVariaveis>::MassMatrix_value(int i,int j)// <== i e j sao modos locais=<
//  {
//    //printf("entrou em MassMatrix_value\n");
//    return((*ptMb_comp)(i,j)*sgn[(i/NumLocalVars)]*sgn[(j/NumLocalVars)])*rho;
//  };
//

// ****************************************************************************
// int PhElem<NumVariaveis>::map(const int & no_local)
// {
//   return (gbnmap[no_local]);
// };

//  // ************************************************************************
//  void PhElem<NumVariaveis>::P_eval_u(const double X[])
//  {
//    int nb=ptr_stdel->nb_val();
//    int Nb=nb*NumLocalVars;
//    int Nn=ptr_stdel->nn_val()*NumLocalVars;
//    int ii,gbi;
//
//    for(int ia=0;ia<NumLocalVars;ia++){
//      for(int i=0;i<nb;i++){
//        ii=i*NumLocalVars+ia;
//        gbi=map(i,ia);
//        u0[ii]=X[gbi]*sgn[i];
//      }
//    }
//    for(int i=Nb;i<Nn;i++){
//      u0[i]=0.0;
//      for(int j=0;j<Nb;j++){
//        u0[i]-=(*ptMcMi_inv)(j,i-Nb)*u0[j];
//      }
//      for(int j=Nb;j<Nn;j++){
//        u0[i]+=(*ptMi_inv)(i-Nb,j-Nb)*b0[j]/rho;
//      }
//    }
//  };
// ******************************************************************
/* void  PhElem::P_eval_u(const double X[],const int & ia)
 {
 int nb=ptr_stdel[ia]->nb_val();
 int nn=ptr_stdel[ia]->nn_val();
 int gbi;
 for(int i=0;i<nb;i++){
 gbi=map(i,ia);
 u0[ia][i]=X[gbi]*sgn[ia][i];
 }
 for(int i=nb;i<nn;i++){
 u0[ia][i]=0.0;
 for(int j=0;j<nb;j++){
 u0[ia][i]-=(*ptMcMi_inv)(j,i-nb)*u0[ia][j];
 }
 for(int j=nb;j<nn;j++){
 u0[ia][i]+=(*ptMi_inv)(i-nb,j-nb)*b0[ia][j]/rho;
 }
 }
 };
 */
// // ******************************************************************
// void PhElem<NumVariaveis>::P_eval_u_p1(const double X[],const int & ia)
// {
//   // Valida para a Equacao_p1 do problema ImiscFlow2F
//   int nb=ptr_stdel[ia]->nb_val();
//   int nn=ptr_stdel[ia]->nn_val();
//   int ii,jj,gbi;
//   for(int i=0;i<nb;i++){
//     gbi=map(i,ia);
//     u0[ia][i]=X[gbi]*sgn[ia][i];
//   }
//   for(int i=nb;i<nn;i++){
//     u0[ia][i]=0.0;
//     for(int j=0;j<nb;j++){
//       u0[ia][i]-=(*ptKcKi_inv)(j,i-nb)*u0[ia][j];
//     }
//     for(int j=nb;j<nn;j++){
//       u0[ia][i]+=(*ptKi_inv)(i-nb,j-nb)*b0[ia][j];
//     }
//   }
// };


// ****************************************************************************
//int PhElem<NumVariaveis>::show_sgn(const int & ia, const int & i)
//{
// return (sgn[ia][i]);
// };

//  double PhElem<NumVariaveis>::show_u(const int & ia,const int & i)
//  {
//    return u0[ia][i];
//  };

// ****************************************************************************
// Calcula os valores da funcao de interpolacao nos vertices
// ****************************************************************************
//  void PhElem<NumVariaveis>::P_eval_vert(const double X[],double f_vert[])
//  {
//    int i,m1;
//    int nb=ptr_stdel->nb_val();
//    int nn=ptr_stdel->nn_val();
//    int ni=nn-nb;
//    int p, q;
//    for(i=0;i<nb;i++)u0[i]=X[gbnmap[i]]*sgn[i];
//    for(i=nb;i<nn;i++){
//      u0[i]=0.0;
//      for(int j=0;j<nb;j++){
//        u0[i]-=(*ptMcMi_inv)(j,i-nb)*u0[j];
//      }
//      for(int j=nb;j<nn;j++){
//        u0[i]+=(*ptMi_inv)(i-nb,j-nb)*b0[j]/rho;
//      }
//    }
//    ptr_stdel->computeVertice(f_vert,u0,ptvert,gbnmap);
//  };

// ***************************************************************************
// void PhElem<NumVariaveis>::P_eval_print(const double X[],const int & numf, FILE *file,
// 			double (*func)(double,double,double))
// {
//   int nn=ptr_stdel[numf]->nn_val();
//   int i;
//   double utemp[nn];
//   P_eval_u(X,numf);
//   for(i=0;i<nn;i++){
//     utemp[i]=u0[numf][i];
//     //printf("P_eval_print var=%d u[%d]=%g\n",numf,i,utemp[i]);
//   }
//   ptr_stdel[numf]->printtofile(file,utemp,func,ptvert,Vert_map);
// };


// ***************************************************************************
//  void PhElem<NumVariaveis>::print_matrices(FILE * fout)
//  {
//    int i,j;
//    int nb=ptr_stdel->nb_val();
//    int nn=ptr_stdel->nn_val();
//
//    // fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mass\n");
//    // for(i=0;i<Nb;i++){
//    //   //fprintf(fout,"linha %d\n",i);
//    //   for(j=0;j<Nb;j++){
//    //     fprintf(fout,"%11.4e ",ptr_stdel->mass(i,j)*rho);
//    //   }
//    //   fprintf(fout,"\n");
//    // }
//
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mb_comp\n");
//    for(i=0;i<nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      for(j=0;j<nb;j++){
//        fprintf(fout,"%11.4e ",(*ptMb_comp)(i,j)*rho);
//      }
//      fprintf(fout,"\n");
//    }
//
//    //  fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mc_transposta\n");
//    //   for(i=0;i<nn-nb;i++){
//    //     //fprintf(fout,"linha %d\n",i);
//    //     for(j=0;j<nb;j++){
//    //       fprintf(fout,"%11.4e ",Mc_T[i][j]);
//    //     }
//    //     fprintf(fout,"\n");
//    //   }
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz Mi_inv\n");
//    for(i=0;i<nn-nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      for(j=0;j<nn-nb;j++){
//        fprintf(fout,"%11.4e ",(*ptMi_inv)(i,j)/rho);
//      }
//      fprintf(fout,"\n");
//    }
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices Matriz McMi_inv\n");
//    for(i=0;i<nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      for(j=0;j<nn-nb;j++){
//        fprintf(fout,"%11.4e ",(*ptMcMi_inv)(i,j));
//      }
//      fprintf(fout,"\n");
//    }
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices vetor b\n");
//    for(i=0;i<nn;i++){
//      //fprintf(fout,"linha %d\n",i);
//      fprintf(fout,"%11.4e ",b0[i]);
//    }
//    fprintf(fout,"\n");
//
//    fprintf(fout,"PhElem<NumVariaveis>::print_matrices gbnmap e sgn  \n");
//    for(i=0;i<nb;i++){
//      //fprintf(fout,"linha %d\n",i);
//      fprintf(fout,"%4d %4d %4d\n",i,gbnmap[i], sgn[i]);
//    }
//    fprintf(fout,"\n");
//  };

// *********************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::compute_JV(const int & ia)
{
  if(ia < NumVariaveis) {
    int q0=ptr_stdel[ia]->NGQP_val();
    //printf("Calculo do Jacobiano: variavel %d\n",ia);
    // Armazenamento dinamico da memoria para JV;
    // Opcao 1: Uma matriz com 3 indices
    //   JV = new double ** [q0];
    //   for(int i=0;i<q0;i++){
    //     JV[i] = new double * [q1];
    //     for(int j=0;j<q1;j++){
    //       JV[i][j] = new double  [q2];
    //     }
    //   }
    
    //JV[ia] = new double [q0*q1*q2]; // opcao 2: um unico vetor
    JV = new double [q0];
    
    //printf("ANTES: Calculo do Jacobiano da var %d  dimensao=%d\n",ia,q0*q1*q2);
    
    ptr_stdel[ia]->Jacobian(ptvert,Vert_map,JV);
  }
  else cout << " elemento nao contem a variavel " << ia << endl;
  //printf("DEPOIS: Calculo do Jacobiano dimensao=%d\n",q0*q1*q2);
};

// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_funcoes_na_borda(EDGE * border)
{
  // Esta procedure permite que se calcule os tracos
  // diretamente das funcoes sem fazer a operacao de
  // fatiamento. Pode ser usada com pontos de Gauss
  // totalmente internos.
  // ***************************************
  // Calcular Phi e seu gradiente na borda *
  // Nao deve usar eval_Phi() pois essa    *
  // usa os pontos de Gauss do elemento,   *
  // que podem nao existir nas bordas      *
  // ***************************************
  int h,i,j,pos;
  int ndim = ptr_stdel[0]->ndim_val();
  double n_e[ndim];
  
  for(i=0;i<NumVariaveis;i++) {//i= variavel
    
    int nborder= ptr_stdel[i]->nborder_val();
    int nn     = ptr_stdel[i]->nn_val();
    int qmax   = ptr_stdel[i]->qborder_val();
    int NGQP   = ptr_stdel[i]->NGQP_val();
    double phi[NGQP];
    
    //	double TP[nborder][nn][qmax];
    
    double *** TP = new double ** [nborder];
    for(h=0;h<nborder; h++) {
      TP[h] = new double * [nn];
      for(j=0;j<nn;j++){
        TP[h][j] = new double [qmax];
      }
    }
    
    //int s0 = nn*ndim*qmax;
    //int s1 =    ndim*qmax;
    
    // No elemento
    for(j=0;j<nn;j++){
      ptr_stdel[i]->eval_Phi(j,phi);  // <-- No elemento
      // Calcular as derivadas de Phi_j
      ptr_stdel[i]->Gradiente(GradPhi[i][j],phi,ptvert,Vert_map);
    }
    
    // *************************************************
    // Calcula os Laplacianos de Phi no elemento padrao
    // apos o calculo do Gradiente
    // *************************************************
    // Aloca memoria para temp
    
    double grad[ndim][NGQP];
    double *ptr_grad[ndim];
    for(int k=0;k<ndim;++k){
      ptr_grad[k]=grad[k];
    }
    // Calcula o Laplaciano
    for(j=0;j<nn;++j){
      for(int l=0;l<NGQP;++l) LaplacianoPhi[i][j][l] = 0.0;
      for(int k=0;k<ndim;k++) { // direcao k
        // Calcular as derivadas de GradPhi[i][j][k]
        ptr_stdel[i]->Gradiente(ptr_grad,GradPhi[i][j][k],ptvert,Vert_map);
        for(int l=0;l<NGQP;l++) LaplacianoPhi[i][j][l] += grad[k][l];
      }
    }
    
    // *********************************************
    // Calcula os tracos de Phi e de seu gradiente
    // no elemento padrao
    // *********************************************
    ptr_stdel[i]->elem_traces(ptvert,Vert_map,sinal,TP,TrGradPhi[i],Jb);
    
    // Implementar nos elementos padroes (Triangle, Linear, LinearLeg, Quadrilateral e Tetrahedral)
    
    for(h=0;h<nborder;h++){
      for(int ndir=0;ndir<ndim;ndir++){
        n_e[ndir]=border[border_num[h]].normal[ndir];
      }
      for(j=0;j<nn;j++){
       	if(ptr_stdel[i]->is_on_border(j,h,pos)){
          for(int q=0;q<qmax;q++){
            TrPhi[i][h][pos][q] = TP[h][j][q];
          }
        }
        
        //	int ini = h*s0 + j*s1;
        for(int q=0;q<qmax;q++){
          // printf("q %3d/%3d ",q,qmax);
          double temp = 0.0;
          for(int ndir=0; ndir < ndim;ndir++){
            //	    TrGradPhi[i][h][j][ndir][q] = TGP[ini+ndir*qmax+q];
            temp+= (perm[ndir]*n_e[ndir])*TrGradPhi[i][h][j][ndir][q]; // gphi_[ndir][q]);
            // printf("Concluido TrGradPhi[%d][%d][%d][%d][%d] = % g\n",i,h,j,ndir,q,TrGradPhi[i][h][j][ndir][q]);
          }
          TrKgradPhi_n[i][h][j][q]=temp;
          //printf("\nTrKgradPhi_n[%d][%d][%d][%d] = % g\n\n",i,h,j,q,TrKgradPhi_n[i][h][j][q]);
        }
      }
    }
    
    // Libera memoria dinamica de TP
    for(h=0;h<nborder; h++) {
      for(j=0;j<nn;j++){
        delete [] TP[h][j]; TP[h][j]=nullptr;
      }
      delete [] TP[h]; TP[h] = nullptr;
    }
    delete [] TP; TP =nullptr;
  }
};

// ***************************************************************************
// ***************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_tracos( EDGE * border)
{
  // OBSOLETA!!!
  // Esta procedure necessita que haja pontos de Gauss nas bordas.
  // Para pontos de Gauss totalmente internos ou casos gerais usar
  // inicia_funcoes_na_borda (ver acima)
  // ******************************
  // Calcular Phi e seu gradiente *
  // ******************************
  int i,j,h,pos;
  int ndim = ptr_stdel[0]->ndim_val();
  // int NGQP = ptr_stdel[0]->NGQP_val();
  double n_e[ndim];
  
  for(i=0;i<NumVariaveis;i++) {//i= variavel
    int nn=ptr_stdel[i]->nn_val();
    int nborder=ptr_stdel[i]->nborder_val();
    int NGQP=ptr_stdel[i]->NGQP_val();
    int qmax = ptr_stdel[i]->qborder_val();
    double phi[NGQP];
    double gphi_[ndim][qmax];
    
    for(j=0; j < nn; j++) {// j=modo
      // calcular os valores de Phi_j nos pontos de quadratura
      ptr_stdel[i]->eval_Phi(j,phi);
      // Calcular as derivadas de Phi_j
      ptr_stdel[i]->Gradiente(GradPhi[i][j],phi,ptvert,Vert_map);
      
      // calcula os tracos de gradiente de Phi
      for(h=0; h < nborder; h++) { // h = aresta
        //printf("Ponto teste em inicia_tracos i=%d j=%d h=%d \n",i,j,h);
        for(int ndir=0; ndir<ndim;ndir++){
          // trace eh válido somente quando ha pontos de Gauss na borda
          ptr_stdel[i]->trace(h,qmax,sinal[h],GradPhi[i][j][ndir],gphi_[ndir]);
          n_e[ndir]=border[border_num[h]].normal[ndir];
        }
        for(int q=0; q < qmax; q++){
          // printf("q %3d/%3d ",q,qmax);
          double temp = 0.0;
          for(int ndir=0; ndir < ndim;ndir++){
            temp += (perm[ndir]*n_e[ndir]*gphi_[ndir][q]);
            TrGradPhi[i][h][j][ndir][q]=gphi_[ndir][q];
          }
          TrKgradPhi_n[i][h][j][q]=temp;
          //printf("Concluido\n");
        }
        // armazena o traco de phi
        if(ptr_stdel[i]->is_on_border(j,h,pos)){
          // printf("armazena traco de phi i=%d j= %d h=%d pos= %d ...\n",i,j,h,pos);
          // trace eh válido somente quando ha pontos de Gauss na borda
          ptr_stdel[i]->trace(h,qmax,sinal[h],phi,TrPhi[i][h][pos]);
          
          //printf(" salvou\n ");
          
        }
        //printf("Concluido h =%d (ne)\n",h);
      }
      //printf(" Concluido j = %d \n",j);
    }
    //printf(" Concluido i = %d \n",i);
  }
};

// ******************************************************************************
// ******************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::set_part_num(const int & num)
{
  if(num < 0) { //calcula
    // ****************************************************
    // Elemento pertence a particao de menor numero
    // dentre as particoes de seus vertices
    // ****************************************************
    int temp;
    part_num=ptvert[Vert_map[0]].part_num;
    for(int i=1;i<numv;i++){
      temp=ptvert[Vert_map[i]].part_num;
      if(part_num > temp) part_num=temp;
    }
  }
  else // part_num eh dado na chamada do metodo
    part_num=num;
};

// ****************************************************************
// Calcula os tracos de funcoes somando os tracos dos modos
// O esforco computacional eh proporcional ao numero de modos
// ****************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::calcula_tracos(Fluids fls)
{
  
  const int ndim=ptr_stdel[0]->ndim_val();
  const int nborder=ptr_stdel[0]->nborder_val();
  const int qborder=ptr_stdel[0]->qborder_val();
  
  const int ns=ptr_stdel[sat]->nn_val();
  const int np=ptr_stdel[pres]->nn_val();
  
  int pos;
  double aux;
  
  // Trsn
  for(int h=0; h < nborder; h++){
    for(int q=0; q < qborder; q++){
      aux=0.0;
      for(int n=0; n < ns; n++){
        if(ptr_stdel[sat]->is_on_border(n,h,pos)){
          aux+= TrPhi[sat][h][pos][q]*u0[sat][n];
        }
      }
      Trsn[h][q]=aux;
    }
  }
  
  // Trpw
  for(int h=0; h < nborder; h++){
    for(int q=0; q < qborder; q++){
      aux=0.0;
      for(int n=0; n < np ;n++){
        if(ptr_stdel[pres]->is_on_border(n,h,pos)){
          aux+= TrPhi[pres][h][pos][q]*u0[pres][n];
        }
      }
      Trpw[h][q]=aux;
    }
  }
  
  // TrKgrad_sn e TrKgrad_pc
  
  for(int h=0; h < nborder; h++){
    for(int q=0; q < qborder; q++){
      double aux1= fls.dpc(Trsn[h][q]);   // <---- fonte de dificuldade
      for(int ndir=0; ndir < ndim; ndir++){
        aux=0.0;
        for(int n=0; n < ns; n++){
          aux+=  TrGradPhi[sat][h][n][ndir][q]*u0[sat][n];
        }
        TrKgrad_sn[h][ndir][q]=aux*perm[ndir];
        TrKgrad_pc[h][ndir][q]=aux1*aux;
      }
    }
  }
  
  // TrKgrad_pw
  
  for(int h=0; h < nborder; h++){
    for(int q=0; q < qborder; q++){
      for(int ndir=0; ndir < ndim; ndir++){
        aux=0.0;
        for(int n=0; n < np; n++){
          aux += TrGradPhi[pres][h][n][ndir][q]*u0[pres][n];
        }
        TrKgrad_pw[h][ndir][q]=aux*perm[ndir];
      }
    }
  }
  
  /*
   //2
   
   // *************************************
   // Passar os tracos de sn e pw para os
   // vetores globais
   // *************************************
   if(gbtrpw!=NULL && gbtrsn!= NULL){
   int ind = gbtrbmap;
   for(int h=0;h<nborder;++h) {
			for(int q=0;q<qborder;++q) {
   gbtrpw[ind]=Trpw[h][q];
   gbtrsn[ind]=Trsn[h][q];
   ind++;
			}
   }
   }
   */
}
// *****************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::teste_gradiente()
{
  const int nn=ptr_stdel[0]->nn_val();
  const int ndim=ptr_stdel[0]->ndim_val();
  // const int nborder=ptr_stdel[0]->nborder_val();
  const int NGQP=ptr_stdel[0]->NGQP_val();
  // const int ns=ptr_stdel[sat]->nn_val();
  // const int np=ptr_stdel[pres]->nn_val();
  double ** der = new  double * [ndim];
  double ** grad = new  double * [ndim];
  for(int k=0;k<ndim;k++) {
    der[k] = new double[NGQP];
    grad[k] = new double[NGQP];
  }
  for(int i=0;i<2;i++){
    for(int j=0; j < nn; j++) {// j=modo
      // calcular os valores de Phi_j nos pontos de quadratura
      double phi[NGQP];
      ptr_stdel[i]->eval_Phi(j,phi);
      // Calcular as derivadas de Phi_j
      ptr_stdel[i]->Gradiente(grad,phi,ptvert,Vert_map);
      ptr_stdel[i]->eval_GradPhi(ptvert,Vert_map,j,der);
      for(int q=0;q<NGQP;q++){
        for(int ndir=0;ndir<ndim;ndir++){
          double aux= grad[ndir][q] - der[ndir][q];
          if(abs(aux) > 1e-14 ){
            printf("var = %2d modo = %2d ponto = %2d dir = %2d  grad = %16.9e der = %16.9e diferenca = %16.9e\n",i,j,q,ndir,grad[ndir][q],der[ndir][q], aux);
          }
        }
      }
    }
  }
  for(int k=0;k<ndim;k++) {
    delete [] der[k]; der[k]=nullptr;
  }
  delete [] der; der=nullptr;
};
// *********************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::echo_traco(FILE * f_eco)
{
  int pos;
  for(int i=0;i<NumVariaveis;++i){
    int nborder=ptr_stdel[i]->nborder_val();
    int nn=     ptr_stdel[i]->nn_val();
    int qmax =  ptr_stdel[i]->qborder_val();
    int ndim =  ptr_stdel[i]->ndim_val();
    
    for(int h=0;h<nborder;h++){
      
      for(int j=0;j<nn;j++){
        
        if(ptr_stdel[i]->is_on_border(j,h,pos)){
          for(int q=0;q<qmax;++q){
            printf("\nTrPhi[%d][%d][%d(%d)][%d] = % g\n\n",i,h,pos,j,q,TrPhi[i][h][pos][q]);
          }
        }
        for(int q=0;q<qmax;++q){
          // printf("q %3d/%3d ",q,qmax);
          // double temp = 0.0;
          for(int ndir=0; ndir < ndim;++ndir){
            printf("Concluido TrGradPhi[%d][%d][%d][%d][%d] = % g\n",i,h,j,ndir,q,TrGradPhi[i][h][j][ndir][q]);
          }
          printf("\nTrKgradPhi_n[%d][%d][%d][%d] = % g\n\n",i,h,j,q,TrKgradPhi_n[i][h][j][q]);
        }
      }
    }
  }
  
  if(f_eco != NULL){
    // *************************************************************
    // Echo dos resultados
    // *************************************************************
    int NGQP=ptr_stdel[0]->NGQP_val();
    fprintf(f_eco,"Jacobiano\n");
    for(int i=0;i<NGQP;++i)fprintf(f_eco,"%11.4e\n",JV[i]);
    
    fprintf(f_eco,"\nPhi para o elemento\n");
    for(int i=0;i<NumVariaveis;++i){
      fprintf(f_eco,"\n variavel %d\n",i);
      int nn     =ptr_stdel[i]->nn_val();
      NGQP       =ptr_stdel[i]->NGQP_val();
      int nborder=ptr_stdel[i]->nborder_val();
      int qmax   =ptr_stdel[i]->qborder_val();
      
      //int ndim =  ptr_stdel[i]->ndim_val();
      
      for(int h=0;h<nborder;++h){
        fprintf(f_eco,"\naresta local %d sinal = %d\n",h,sinal[h]);
        for(int j=0;j<nn;++j){
          
          if(ptr_stdel[i]->is_on_border(j,h,pos)){
            fprintf(f_eco,"modo %d: ",j);
            for(int q=0;q<qmax;++q){
              //fprintf(f_eco,"\nTrPhi[%d][%d][%d(%d)][%d] = % g\n\n",i,h,pos,j,q,TrPhi[i][h][pos][q]);
              fprintf(f_eco,"%g ",TrPhi[i][h][pos][q]);
            }
            fprintf(f_eco,"\n");
          }
          /*
           for(int q=0;q<qmax;++q){
           // printf("q %3d/%3d ",q,qmax);
           // double temp = 0.0;
           for(int ndir=0; ndir < ndim;++ndir){
           fprintf(f_eco,"Concluido TrGradPhi[%d][%d][%d][%d][%d] = % g\n",i,h,j,ndir,q,TrGradPhi[i][h][j][ndir][q]);
           }
           fprintf(f_eco,"\nTrKgradPhi_n[%d][%d][%d][%d] = % g\n\n",i,h,j,q,TrKgradPhi_n[i][h][j][q]);
           }
           */
        }
      }
      /*
       for(int j=0;j<nn;++j){
       fprintf(f_eco,"         modo %d\n",j);
       fprintf(f_eco,"Grad_Phi[%d][%d]\n",i,j);
       for(int m=0;m<NGQP;++m){
       fprintf(f_eco," m= %3d %11.4e  %11.4e\n",m,GradPhi[i][j][0][m],GradPhi[i][j][1][m]);
       }
       //printf("AQUI1 ***************************************************\n");
       for(int h=0;h<nborder;++h){
       fprintf(f_eco,"  Traco na aresta %d \n",h);
       for(int m=0;m<qmax;++m){
       fprintf(f_eco,"%11.4e %11.4e\n",TrGradPhi[i][h][j][0][m],TrGradPhi[i][h][j][1][m]);
       }
       }
       //printf("AQUI2 **************************************************\n");
       }
       */
    }
    // **********************************************************
    // Fim do echo dos tracos de Phi e GradPhi
    // **********************************************************
  }
  
};
/*
 template<int NumVariaveis>
 void PhElem::construir_PhiArray()
 {
 int NGQP = ptr_stdel[0]->NGQP_val();
 int nn =  ptr_stdel[0]->nn_val();
 
 double phi[NGQP];
 PhiArray = new double [nn*NGQP];
 
 int count = 0;
 for(int j=0; j < nn; j++) {// j=modo
 // calcular os valores de Phi_j nos pontos de quadratura
 ptr_stdel[0]->eval_Phi(j,phi);
 for(int i=0;i<NGQP;i++) PhiArray[count++]=phi[i];
 }
 };
 template<int NumVariaveis>
 const double * PhElem::eval_Phi(const int i)
 {
 int NGQP = ptr_stdel[0]->NGQP_val();
 return (PhiArray + (i*NGQP));
 };
 */

template<int NumVariaveis>
void PhElem<NumVariaveis>::vetor_superficie(const int & num_local,double & area, double normal[3])
{
  ptr_stdel[0]->superficie_externa(Vert_map,ptvert,num_local,area,normal);
};
// ****************************************************************************
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::perm_val(double K[3])
{
  for(int i=0;i<3;i++)K[i]=perm[i];
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_sn(const int & h,double * saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  for(int q=0;q<qmax;q++)saida[q]=Trsn[h][q];
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_pw(const int & h,double *saida)
{
  int qmax=ptr_stdel[1]->qborder_val();
  for(int q=0;q<qmax;q++)saida[q]=Trpw[h][q];
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_Kgrad_pw(const int & h,double ** saida)
{
  int qmax=ptr_stdel[1]->qborder_val();
  int ndim = ptr_stdel[1]->ndim_val();
  for(int k=0;k<ndim;k++){
    for(int q=0;q < qmax;q++){
      saida[k][q]=TrKgrad_pw[h][k][q];
    }
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_Kgrad_pc(const int & h,double ** saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  int ndim = ptr_stdel[1]->ndim_val();
  for(int k=0;k<ndim;k++){
    for(int q=0;q < qmax;q++) {
      saida[k][q]=TrKgrad_pc[h][k][q];
    }
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_Kgrad_sn(const int & h,double ** saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  int ndim = ptr_stdel[1]->ndim_val();
  
  for(int k=0;k < ndim;k++){
    for(int q=0;q < qmax;q++){
      if(std::isnan(TrKgrad_sn[h][k][q]))cout << "Float was Not a Number: TrKgrad " << q << endl;
      saida[k][q]=TrKgrad_sn[h][k][q];
    }
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_phi(const int & lado,const int & ivar,
                                           const int & ind, double * saida)
{
  int qmax=ptr_stdel[ivar]->qborder_val();
  //int NGQP=ptr_stdel[ivar]->NGQP_val();
  //double phi[NGQP];
  int pos;
  //ptr_stdel[ivar]->eval_Phi(ind,phi);
  //ptr_stdel[ivar]->trace(lado,qmax,sinal[lado],phi,saida);
  if(ptr_stdel[ivar]->is_on_border(ind,lado,pos)){
    for(int q=0;q<qmax;q++){
      saida[q]=TrPhi[ivar][lado][pos][q];
    }
  }
  else {
    for(int q=0;q<qmax;q++){
      saida[q]=0.0;
    }
  }
};
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_grad_phi(const int & lado,const int & ivar,
                                                const int & ind,double ** saida)
{// Ficou desnecessaria quando criei Traco_Kgrad_phi_n abaixo
  int qmax=ptr_stdel[0]->qborder_val();
  int ndim = ptr_stdel[1]->ndim_val();
  for(int k=0;k<ndim;k++){
    for(int q=0;q<qmax;q++)
      saida[k][q]=TrGradPhi[ivar][lado][ind][k][q];
  }
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::Traco_Kgrad_phi_n(const int & lado,const int & ivar,
                                                   const int & ind,double * saida)
{
  int qmax=ptr_stdel[0]->qborder_val();
  for(int q=0;q<qmax;q++)
    saida[q]=TrKgradPhi_n[ivar][lado][ind][q];
};
// ****************************************************************************
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::VolumeIntegrals_IG(Fluids fls,  int & count,
                                                    int * Ti,int * Tj, double * Tx,
                                                    double * B)
{
  //cout << "VolumeIntegrals_IG  \n";
  
  int qmax=ptr_stdel[sat]->qborder_val();
  int NGQP=ptr_stdel[sat]->NGQP_val();
  int ndim=ptr_stdel[sat]->ndim_val();
  int nborder  =ptr_stdel[sat]->nborder_val();
  
  int np;//ns
  int h,k,m;
  double aux;
  double res0[NGQP];
  double res[ndim][NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  for(k=0;k<ndim;k++){
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP){
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  for(m=0;m<NGQP;m++)JW[m]*=JV[m];
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++)pc[m]=fls.pressao_capilar(sn[m]);
  ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
  
  for(m=0;m<NGQP;m++){
    // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
    for(k=0;k<ndim;k++){// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
    }
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
  }
  
  
  //printf("VolumeIntegrals_IG Ate aqui 1 ...\n");
  
  /*
   // salva os tracos de K * gradientes
   calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
   double test[nborder][ndim][qmax];
   for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   test[h][k][q]=TrKgrad_sn[h][k][q];
   }
   }
   }
   */
  
  for(h=0;h<nborder;h++){
    ptr_stdel[sat] ->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++){
      ptr_stdel[sat] ->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat] ->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  
  
  /*
   for(h=0;h<nborder;h++)
   {
   for(k=0;k<ndim;k++)
   {
   for(int q=0;q<qmax;q++)
   {
   if(test[h][k][q]!=TrKgrad_pw[h][k][q])
   {
   printf("Em h= %d ndir =%d q=%d test = %g diferenca de TrKgrad_sn = %g   \n", h,k,q,test[h][k][q],
   TrKgrad_sn[h][k][q]-test[h][k][q]);
   }
   }
   }
   }
   */
  
  //
  // printf("VolumeIntegrals_IG Ate aqui 2\n");
  np=ptr_stdel[pres]->nn_val();
  // ns=ptr_stdel[sat ]->nn_val();
  for(int rp=0;rp<np;rp++){ // PE_VI
    // *****************************
    // calculo do vetor
    // *****************************
    
    // integral 2 com sinal trocado
    for(int k=0; k < ndim; k++) {
      prodvv(NGQP,Kgrad_pc[k],GradPhi[pres][rp][k],res[k]);
    }
    for(int k=1; k < ndim; k++) {
      somavv(NGQP,res[0],res[k],res[0]);
    }
    prodvv(NGQP,lambdan,res[0],res0);
    integral(NGQP,JW,res0,aux);
    if(std::isnan(aux))cout << "VolumeIntegrals_IG  Float was Not a Number: aux " << aux << endl;
    B[gbnmap[pres][rp]] -= aux; // primeiro valor de B
    
    // Fim do vetor em PE_VI
    
    // *****************************
    // Matriz
    // *****************************
    for(int lp=0;lp<np;lp++){ // PE_VI_PP
      for(int k=0; k < ndim; k++) {
        prodev(NGQP,perm[k],GradPhi[pres][lp][k],res[k]);
      }
      for(int k=0; k < ndim; k++) {
        prodvv(NGQP,res[k],GradPhi[pres][rp][k],res[k]);
      }
      for(int k = 1; k < ndim; k++) {
        somavv(NGQP,res[0],res[k],res[0]);
      }
      prodvv(NGQP,lambdat,res[0],res0);
      integral(NGQP,JW,res0,aux);
      if(std::isnan(aux))cout << "VolumeIntegrals_IG  Float was Not a Number 2: aux " << aux << endl;
      Ti[count]=gbnmap[pres][rp];
      Tj[count]=gbnmap[pres][lp];
      Tx[count]=aux;
      count++;
    }
  }
  // OK 22/05/2008
};

//   // ****************************************************************************
//  // ****************************************************************************
// template<int NumVariaveis>
//  void PhElem<NumVariaveis>::VolumeIntegrals_IG_Epetra(Fluids fls,
//  		       Epetra_FECrsMatrix & A,
//  		       Epetra_FEVector & RHS)
//  {
//   cout << "PhElem<NumVariaveis>::VolumeIntegrals_IGa \n";
//  #define indice(n0,var,i) ((n0*var)+i)
//  #define npos(ntot,i,j) ((ntot*i)+j)
//
//
//    int qmax=ptr_stdel[sat]->Q_val(0);
//    int NGQP=ptr_stdel[sat]->Q_val(0)*ptr_stdel[sat]->Q_val(1);
//
//    int ns,np;
//    int h,i,j,k,m,ivar;
//    double aux,aaux;
//    double res0[NGQP],res1[NGQP],res2[NGQP],res3[NGQP];
//    double JW[NGQP];
//    double sn[NGQP], pw[NGQP],pc[NGQP];
//    double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
//    double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
//    double d_pc[NGQP],d2_pc[NGQP];
//    double phi_r[NGQP],phi_l[NGQP];
//    double grad_sn[3][NGQP],grad_pw[3][NGQP],grad_pc[3][NGQP];
//    double *ptr_gsn[3],*ptr_gpw[3],*ptr_gpc[3];
//    double mun=fls.show_mun();
//    double muw=fls.show_muw();
//    for(k=0;k<3;k++){
//      ptr_gsn[k]=grad_sn[k];
//      ptr_gpw[k]=grad_pw[k];
//      ptr_gpc[k]=grad_pc[k];
//    }
//
//    int ndim=ptr_stdel[sat]->ndim_val();
//    int ne  =ptr_stdel[sat]->ne_val();
//
//    // calcula os pesos de Gauss multiplicados pelos jacobianos
//    if(ptr_stdel[sat]->w_vec(JW) != NGQP){
//      printf("Erro de dimensionamento de JW\n");
//      exit(0);
//    }
//    for(m=0;m<NGQP;m++)JW[m]*=JV[m];
//
//    // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
//    ptr_stdel[sat]->evalGQ(sn,u0[sat]);
//    ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
//    ptr_stdel[pres]->evalGQ(pw,u0[pres]);
//    ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
//
//    // Calcula pc e seu gradiente nos pontos de Gauss
//    for(m=0;m<NGQP;m++)pc[m]=fls.pressao_capilar(sn[m]);
//    ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
//
//    for(m=0;m<NGQP;m++){
//      // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
//      for(k=0;k<ndim;k++){// caso especial de tensor diagonal
//        grad_sn[k][m]*=perm[k];
//        grad_pw[k][m]*=perm[k];
//        grad_pc[k][m]*=perm[k];
//      }
//      // Calculo dos escalares
//      aux=sn[m];
//      lambdan[m]=fls.Krn(aux)/mun;
//      lambdaw[m]=fls.Krw(aux)/muw;
//      lambdat[m]=lambdan[m]+lambdaw[m];
//      d_lambdan[m]=fls.dkrn(aux)/mun;
//      d_lambdaw[m]=fls.dkrw(aux)/muw;
//      d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
//      d_pc[m]=fls.dpc(aux);
//      d2_pc[m]=fls.d2pc(aux);
//    }
//
//    // salva os tracos de K * gradientes
//    /*
//    for(h=0;h<ne;h++){
//      ptr_stdel[sat]->trace(h,qmax,sinal[h],sn,Trsn[h]);
//      ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
//      for(k=0;k<ndim;k++){
//        ptr_stdel[sat]->trace(h,qmax,sinal[h],grad_sn[k],TrKgrad_sn[h][k]);
//        ptr_stdel[sat]->trace(h,qmax,sinal[h],grad_pc[k],TrKgrad_pc[h][k]);
//        ptr_stdel[pres]->trace(h,qmax,sinal[h],grad_pw[k],TrKgrad_pw[h][k]);
//      }
//    }
//    */
//    calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
//    np=ptr_stdel[pres]->nn_val();
//    ns=ptr_stdel[sat ]->nn_val();
//    for(int rp=0;rp<np;rp++){ // PE_VI
//      // *****************************
//      // calculo do vetor
//      // *****************************
//
//      // integral 2 com sinal trocado
//      prodgg(NGQP,grad_pc[0],grad_pc[1],
//  	   GradPhi[pres][rp][0],GradPhi[pres][rp][1],res0);
//      prodvv(NGQP,lambdan,res0,res0);
//      integral(NGQP,JW,res0,aux);
//      //B[gbnmap[pres][rp]] -= aux; // primeiro valor de B
//      B[rp] -= aux;
//
//      // Fim do vetor em PE_VI
//
//      // *****************************
//      // Matriz
//      // *****************************
//      for(int lp=0;lp<np;lp++){ // PE_VI_PP
//        prodev(NGQP,perm[0],GradPhi[pres][lp][0],res0);
//        prodev(NGQP,perm[1],GradPhi[pres][lp][1],res1);
//        prodgg(NGQP,res0,res1,GradPhi[pres][rp][0],GradPhi[pres][rp][1],res0);
//        prodvv(NGQP,lambdat,res0,res0);
//        integral(NGQP,JW,res0,aux);
//        //Ti[count]=gbnmap[pres][rp];
//        //Tj[count]=gbnmap[pres][lp];
//        mx[ ( (rp * np) + lp) ] = aux;
//        count++;
//      }
//    }
//    for(int rp=0;rp<np;rp++){
//      indx[rp] = gbnmap[pres][rp];
//    }
//
//    A.InsertGlobalValues(ntot,indx,mx,Epetra_FECrsMatrix::ROW_MAJOR);
//    RHS.SumIntoGlobalValues(ntot,indx,B);
//    delete [] indx;indx=0;
//    delete [] mx; mx=0;
//    delete [] B; B=0;
//    // OK 22/05/2008
//  };

// ****************************************************************************
// Volume Integrals for regular elements creating UMFPACK arrays
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::VolumeIntegrals_UMFPACK(const double Dt,Fluids fls, int & count,
                                                         int * Ti,int * Tj, double * Tx,
                                                         double * B,
                                                         double * gbtrsn, double * gbtrpw)
{
  //cout << "DG_VI.cc PhElem<NumVariaveis>::VolumeIntegrals_UMFPACK\n";
  
  const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  const int ndim=ptr_stdel[sat]->ndim_val();
  const int nborder=ptr_stdel[sat]->nborder_val();
  
  int h,k,m;
  double aux,aaux;
  double res0[NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  double phi_r[NGQP],phi_l[NGQP]; // em PhElem
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  
  
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  
  for(k=0;k<ndim;k++) {
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  // calcula os pesos de Gauss no vetor JW
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  
  // multiplica os pesos de Gauss pelo Jacobiano JV
  for(m=0;m<NGQP;m++){
    JW[m]*=JV[m];
    // printf(" JW  = %g    JV = %g \n", JW[m], JV[m]);
  }
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++) {
    pc[m]=fls.pressao_capilar(sn[m]);
  }
  ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
  
  // ***********************************************************************
  // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
  for(m=0;m<NGQP;m++) {
    for(k=0;k<ndim;k++) {// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
    }
    
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
  }
  
  // *********************************************************************
  // salva os tracos de K * gradientes
  
  
  /*
   //calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
   double test[nborder][ndim][qmax];
   for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   test[h][k][q]=TrKgrad_pw[h][k][q];
   }
   }
   }
   */
  for(h=0;h<nborder;h++) {
    ptr_stdel[sat]->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++) {
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  /* for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   if(test[h][k][q]!=TrKgrad_pw[h][k][q]) {
   printf("Em h= %d ndir =%d q=%d test = %g diferenca de TrKgrad_sn = %g   \n", h,k,q,test[h][k][q],TrKgrad_sn[h][k][q]-test[h][k][q]);
   }
   }
   }
   } */
  // *************************************
  // Passar os tracos de sn e pw para os
  // vetores globais
  // *************************************
  int ind = gbtrbmap;
  for(h=0;h<nborder;h++) {
    for(int q=0;q<qmax;q++) {
      gbtrpw[ind]=Trpw[h][q];
      gbtrsn[ind]=Trsn[h][q];
      ind++;
    }
  }
  
  const int ns=ptr_stdel[sat ]->nn_val();
  const int np=ptr_stdel[pres]->nn_val();
  const int ntot=ns+np;// Num total de modos no elemento
  const int n0 = ptr_stdel[0]->nn_val();// Num de modos para a variavel 0
  
  // vetor com indices
  int mr,ml;
  double  mx [ntot*ntot];
  for(int i=0;i<ntot*ntot;i++) {
    mx[i]=0.0;
  }
  
  for(int rp=0;rp<np;rp++) { // PE_VI
    mr=indice(n0,pres,rp);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    aux=0.0;
    double btemp;
    for(int q=0;q<NGQP;q++) { // PE_VI
      
      double temp_pw = 0.0;
      double temp_pc = 0.0;
      for (int k = 0; k < ndim; k++) {
        temp_pw += Kgrad_pw[k][q]*GradPhi[pres][rp][k][q];
        temp_pc += Kgrad_pc[k][q]*GradPhi[pres][rp][k][q];
      }
      
      aux+=JW[q]*(lambdat[q]*temp_pw + lambdan[q]*temp_pc);
    }
    
    btemp=aux; // inicializa valor de B em PE_VI
    
    // termo de fonte
    aaux=qn+qw;
    if(aaux!=0.0) {
      ptr_stdel[pres]->eval_Phi(rp,phi_r);
      aux=0.0;
      for(int q=0;q<NGQP;q++)
        aux+=JW[q]*phi_r[q];
      //B[gbnmap[pres][rp]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    B[gbnmap[pres][rp]]+=btemp;
    // Fim do vetor em PE_VI
    
    // *****************************
    // Matriz
    // *****************************
    for(int lp=0;lp<np;lp++) { // PE_VI_PP
      ml=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double temp = 0.0;
        for (int k = 0; k < ndim; k++) {
          temp += perm[k]*GradPhi[pres][lp][k][m]*GradPhi[pres][rp][k][m];
        }
        
        aux+=( temp * lambdat[m]*JW[m]);
      }
      mx[npos(ntot,mr,ml)]+=aux;
    } // Fim de PE_VI_PP
    
    for(int ls=0;ls<ns;ls++) { // PE_VI_SP
      ml=indice(n0,sat,ls);
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double Vaux[ndim];
        for(k=0;k<ndim;k++) {
          Vaux[k]=( d_lambdat[m]*Kgrad_pw[k][m]+
                   d_lambdan[m]*Kgrad_pc[k][m]+
                   lambdan[m]*d2_pc[m]*Kgrad_sn[k][m] ) *
          phi_l[m]+
          lambdan[m]*d_pc[m]*perm[k]*GradPhi[sat][ls][k][m];
        }
        
        aaux=0.0;
        for(k=0;k<ndim;k++) {
          aaux+=Vaux[k]*GradPhi[pres][rp][k][m];
        }
        
        aux+=(aaux*JW[m]);
      }
      mx[npos(ntot,mr,ml)]+=aux;
    }// Fim de PE_VI_SP
  } // Fim de PE_VI
  
  for(int rs=0;rs<ns;rs++) { // SE_VI
    mr=indice(n0,sat,rs);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    double btemp;
    ptr_stdel[sat]->eval_Phi(rs,phi_r);
    
    aux=0.0;
    for(m=0;m<NGQP;m++) {
      aux+=JW[m]*(sn[m]-sna[m])*phi_r[m];
    }
    aux*=(-porosidade/Dt);
    
    // segunda integral
    aaux=0.0;
    for(m=0;m<NGQP;m++) {
      
      double temp = 0.0;
      for (int k =0; k <ndim; k++) {
        temp += Kgrad_pw[k][m]*GradPhi[sat][rs][k][m];
      }
      aaux+=JW[m]*lambdaw[m]*temp;
    }
    
    aux+=aaux;
    //B[gbnmap[sat][rs]]+=aux; // inicializa valor de B em SE_VI
    btemp=aux;
    // termo de fonte
    
    aaux=qw;
    
    if(aaux!=0.0) {
      integral(NGQP,JW,phi_r,aux);
      //B[gbnmap[sat][rs]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    B[gbnmap[sat][rs]] += btemp ;
    // *************************
    // Matriz
    // *************************
    
    for(int lp=0;lp<np;lp++) { // SE_VI_PS
      ml=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double temp = 0.0;
        for (int k = 0; k < ndim; k++) {
          temp += perm[k]*GradPhi[pres][lp][k][m]*GradPhi[sat][rs][k][m];
        }
        
        aux+=(temp * lambdaw[m]*JW[m]);
      }
      mx[npos(ntot,mr,ml)]+=aux;
    }
    
    for(int ls=0;ls<ns;ls++) { // SE_VI_SS
      ml=indice(n0,sat,ls);
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        aux+=JW[m]*(d_lambdaw[m]*res0[m]-porosidade/Dt*phi_r[m])*phi_l[m];
      }
      mx[npos(ntot,mr,ml)]+=aux;
    }
  }
  // ***************************************************************************
  // Passar a matriz mx para os vetores Ti, Tj e Tx
  // ***************************************************************************
  for(int rp=0;rp<np;rp++) {
    mr=indice(n0,pres,rp);
    
    for(int lp=0;lp<np;lp++) {
      ml=indice(n0,pres,lp);
      Ti[count]=gbnmap[pres][rp];
      Tj[count]=gbnmap[pres][lp];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
    
    for(int ls=0;ls<ns;ls++) {
      ml=indice(n0,sat,ls);
      Ti[count]=gbnmap[pres][rp];
      Tj[count]=gbnmap[sat][ls];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
  }
  for(int rs=0;rs<ns;rs++) {
    mr=indice(n0,sat,rs);
    
    for(int lp=0;lp<np;lp++) {
      ml=indice(n0,pres,lp);
      Ti[count]=gbnmap[sat][rs];
      Tj[count]=gbnmap[pres][lp];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
    
    for(int ls=0;ls<ns;ls++) {
      ml=indice(n0,sat,ls);
      Ti[count]=gbnmap[sat][rs];
      Tj[count]=gbnmap[sat][ls];
      Tx[count]= mx[npos(ntot,mr,ml)];
      count++;
    }
  }
  // delete [] mx; mx=nullptr;
};
// alterado em 21/09/2011
// *****************************************************************************

// *****************************************************************************
// Usar se precisar atualizar so o termo de DT nas integrais de volume
//	 da equacao de sn
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::VolumeIntegralsT(const double Dt_new, const double Dt_old,
                                                  int & count,
                                                  int * Ti,int * Tj, double * Tx,
                                                  double * B)
{
  //cout << "DG_VI.cc  PhElem<NumVariaveis>::VolumeIntegralsT\n";
  
  // const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  
  int m;
  double aux;
  double JW[NGQP];
  double sn[NGQP];
  double phi_r[NGQP]; // em PhElem
  
  const double fator=(1/Dt_new - 1/Dt_old);  // Para simples correcao do termo contendo Dt
  
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  
  for(m=0;m<NGQP;m++)JW[m]*=JV[m];
  
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  
  //const int np=ptr_stdel[pres]->nn_val();
  const int ns=ptr_stdel[sat ]->nn_val();
  
  for(int rs=0;rs<ns;rs++) { // SE_VI
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    
    ptr_stdel[sat]->eval_Phi(rs,phi_r);
    
    aux=0.0;
    for(m=0;m<NGQP;m++) {
      aux+=JW[m]*(sn[m]-sna[m])*phi_r[m];
    }
    aux*=(-porosidade*fator);
    
    B[gbnmap[sat][rs]]+=aux; // corrige valor de B em SE_VI
    
    
    // *************************
    // Matriz
    // *************************
    
    for(int ls=0;ls<ns;ls++) { // SE_VI_SS
      Ti[count]=gbnmap[sat][rs];
      Tj[count]=gbnmap[sat][ls];
      Tx[count]=(-porosidade*fator)*Mass_sn[rs][ls];//aux;
      count++;
    }
  }
};

// ****************************************************************************
// Volume Integrals for ghost elements - calcula os tracos
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::VolumeTracos(const double Dt,Fluids fls, double * gbtrsn, double * gbtrpw)
{
  //cout << "DG_VI.cc PhElem<NumVariaveis>::VolumeTacos\n";
  
  const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  const int nborder = ptr_stdel[sat]->nborder_val();
  const int ndim=ptr_stdel[sat]->ndim_val();
  
  int h,k,m;
  double aux;
  //  double res0[NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  // double phi_r[NGQP],phi_l[NGQP];
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  
  for(k=0;k<ndim;k++) {
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  for(m=0;m<NGQP;m++) JW[m]*=JV[m];
  
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++)pc[m]=fls.pressao_capilar(sn[m]);
  ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);
  // ***********************************************************************
  // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
  for(m=0;m<NGQP;m++) {
    for(k=0;k<ndim;k++) {// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
    }
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
  }
  // *********************************************************************
  // salva os tracos de K * gradientes
  //calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
  
  for(h=0;h<nborder;h++) {
    ptr_stdel[sat]->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++) {
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  
  //2
  
  // *************************************
  // Passar os tracos de sn e pw para os
  // vetores globais
  // *************************************
  int ind = gbtrbmap;
  for(h=0;h<nborder;h++) {
    for(int q=0;q<qmax;q++) {
      gbtrpw[ind]=Trpw[h][q];
      gbtrsn[ind]=Trsn[h][q];
      ind++;
    }
  }
};

/******************************************************************************/
/*  Uso de Trilinos      */
/******************************************************************************/
// ****************************************************************************
// Volume Integrals for regular elements creating Epetra_FECrsMatrix and
// Epetra_FEVector
// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::VolumeIntegrals(const double Dt,Fluids fls,
                                                 Teuchos::RCP<Epetra_FECrsMatrix>  A,
                                                 Teuchos::RCP<Epetra_FEVector> RHS,
                                                 double * gbtrsn, double * gbtrpw)
{
  //cout<< "DG_VI.cc  PhElem<NumVariaveis>::VolumeIntegrals\n";
  //#define indice(n0,var,i) ((n0*var)+i)
  //#define npos(ntot,i,j) ((ntot*i)+j)
  
  
  const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  const int nborder=ptr_stdel[sat]->nborder_val();
  const int ndim=ptr_stdel[sat]->ndim_val();
  
  int h,k,m;
  double aux,aaux;
  double res0[NGQP];
  double JW[NGQP];
  double sn[NGQP], pw[NGQP],pc[NGQP];
  double lambdaw[NGQP],lambdan[NGQP],lambdat[NGQP];
  double d_lambdaw[NGQP],d_lambdan[NGQP],d_lambdat[NGQP];
  double d_pc[NGQP],d2_pc[NGQP];
  double phi_r[NGQP],phi_l[NGQP];
  double Kgrad_sn[ndim][NGQP],Kgrad_pw[ndim][NGQP],Kgrad_pc[ndim][NGQP];
  double *ptr_gsn[ndim],*ptr_gpw[ndim],*ptr_gpc[ndim];
  
  for(k=0;k<ndim;k++) {
    ptr_gsn[k]=Kgrad_sn[k];
    ptr_gpw[k]=Kgrad_pw[k];
    ptr_gpc[k]=Kgrad_pc[k];
  }
  
  double mun=fls.show_mun();
  double muw=fls.show_muw();
  // calcula os pesos de Gauss multiplicados pelos jacobianos
  if(ptr_stdel[sat]->w_vec(JW) != NGQP) {
    printf("Erro de dimensionamento de JW\n");
    exit(0);
  }
  for(m=0;m<NGQP;m++)
    JW[m]*=JV[m];
  
  // Calcula as variaveis sn e pw e seus gradientes nos pontos de Gauss
  ptr_stdel[sat]->evalGQ(sn,u0[sat]);
  ptr_stdel[sat]->Gradiente(ptr_gsn,sn,ptvert,Vert_map);
  ptr_stdel[pres]->evalGQ(pw,u0[pres]);
  ptr_stdel[pres]->Gradiente(ptr_gpw,pw,ptvert,Vert_map);
  
  // Calcula pc e seu gradiente nos pontos de Gauss
  for(m=0;m<NGQP;m++){
    pc[m]=fls.pressao_capilar(sn[m]);
  }
  ptr_stdel[sat]->Gradiente(ptr_gpc,pc,ptvert,Vert_map);// <-- derivada por colocacao; ptr_gpc esta no espaco aproximante
  for(m=0;m<NGQP;m++) {
    // Calculo dos escalares
    aux=sn[m];
    lambdan[m]=fls.Krn(aux)/mun;
    lambdaw[m]=fls.Krw(aux)/muw;
    lambdat[m]=lambdan[m]+lambdaw[m];
    d_lambdan[m]=fls.dkrn(aux)/mun;
    d_lambdaw[m]=fls.dkrw(aux)/muw;
    d_lambdat[m]=d_lambdan[m]+d_lambdaw[m];
    d_pc[m]=fls.dpc(aux);
    d2_pc[m]=fls.d2pc(aux);
    for(k=0;k<ndim;k++) {// caso especial de tensor diagonal
      Kgrad_sn[k][m]*=perm[k];
      Kgrad_pw[k][m]*=perm[k];
      Kgrad_pc[k][m]*=perm[k];
      // ***********************************************************************
      // K * Grad de sn, pw e pc; grad_'s ja estao multiplicados por perm
      // ***********************************************************************
      
      // ---------------------------------------------------------------------------------------------------------
      // Essa opcao gera problemas de convergencia
      //Kgrad_pc[k][m]=d_pc[m]*Kgrad_sn[k][m]; // <----- 06/01/2014 \Nabla{pc} = dpc * \Nabla{sn}; nao esta no espaco aproximante
    }
  }
  // *********************************************************************
  // salva os tracos de K * gradientes
  // Usar uma das duas alternativas abaixo
  
  // *******************************************************************
  // Alternativa 1: Calcula os tracos somando os tracos das
  // funcoes elementares previamente calculados
  // O esforco computacional eh proporcional ao numero de modos e ao
  // numero de pontos de Gauss
  //calcula_tracos(fls);// <--- tem problemas no calculo de TrKgrad_pc (usa grad_pc= dpc * grad_sn)
  // *******************************************************************
  /*
   double test[nborder][ndim][qmax];
   for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   test[h][k][q]=TrKgrad_pc[h][k][q];
   }
   }
   }
   */
  // ou
  
  // *******************************************************************
  // Alternativa 2:
  // Interpola as funcoes nos pontos de Gauss nas bordas
  // O esforco computacional eh proporcional ao numero de pontos de Gauss
  // e nao depende do numero de modos
  // *******************************************************************
  
  for(h=0;h<nborder;h++) {
    ptr_stdel[sat] ->trace(h,qmax,sinal[h],sn,Trsn[h]);
    ptr_stdel[pres]->trace(h,qmax,sinal[h],pw,Trpw[h]);
    for(k=0;k<ndim;k++) {
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_sn[k],TrKgrad_sn[h][k]);
      ptr_stdel[sat ]->trace(h,qmax,sinal[h],Kgrad_pc[k],TrKgrad_pc[h][k]);  // <-- tem que ser calculado assim
      ptr_stdel[pres]->trace(h,qmax,sinal[h],Kgrad_pw[k],TrKgrad_pw[h][k]);
    }
  }
  
  // Teste
  /*
   for(h=0;h<nborder;h++) {
   for(k=0;k<ndim;k++) {
   for(int q=0;q<qmax;q++) {
   if(test[h][k][q]!=TrKgrad_pc[h][k][q]) {
   printf("h= %1d q=%1d %11.4e %11.4e %11.4e   \n", h,q,test[h][k][q],TrKgrad_pc[h][k][q],TrKgrad_pc[h][k][q]-test[h][k][q]);
   }
   }
   }
   }
   */
  // Fim do calculo dos tracos
  // ****************************************************************************
  // 3
  
  // *************************************
  // Passar os tracos de sn e pw para os
  // vetores globais
  // *************************************
  // cout << "gbtrbmap = "<< gbtrbmap << "\n";
  int ind = gbtrbmap;
  for(h=0;h<nborder;h++) {
    for(int q=0;q<qmax;q++) {
      gbtrpw[ind]=Trpw[h][q];
      gbtrsn[ind]=Trsn[h][q];
      ind++;
    }
  }
  
  const int ns=ptr_stdel[sat ]->nn_val();
  const int np=ptr_stdel[pres]->nn_val();
  const int ntot=ns+np;// Num total de modos no elemento
  const int n0 = ptr_stdel[0]->nn_val();// Num de modos para a variavel 0
  
  int mi,mj;
  double mx /*= new double*/ [ntot*ntot];
  double B  /*= new double*/ [ntot];
  int  indx /*= new int*/ [ntot];
  for(int i=0;i<ntot*ntot;i++) {
    mx[i]=0.0;
  }
  for(int i=0;i<ntot;i++) {
    B[i]=0.0;
  }
  
  for(int rp=0;rp<np;rp++) { // PE_VI
    mi=indice(n0,pres,rp);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    aux=0.0;
    double btemp;
    for(int q=0;q<NGQP;q++) { // PE_VI
      
      double temp_pw = 0.0;
      double temp_pc = 0.0;
      
      for (int k = 0; k < ndim; k++) {
        temp_pw += Kgrad_pw[k][q]*GradPhi[pres][rp][k][q];
        temp_pc += Kgrad_pc[k][q]*GradPhi[pres][rp][k][q];
      }
      aux+=JW[q]*(lambdat[q]*temp_pw + lambdan[q]*temp_pc);
    }
    
    btemp=aux; // inicializa valor de B em PE_VI
    
    // termo de fonte
    aaux=qn+qw;
    if(aaux!=0.0) {
      ptr_stdel[pres]->eval_Phi(rp,phi_r);
      aux=0.0;
      for(int q=0;q<NGQP;q++)
        aux+=JW[q]*phi_r[q];
      //B[gbnmap[pres][rp]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    //B[gbnmap[pres][rp]]+=btemp;
    B[mi] += btemp;
    // Fim do vetor em PE_VI
    
    // *****************************
    // Matriz
    // *****************************
    for(int lp=0;lp<np;lp++) { // PE_VI_PP
      mj=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        double temp = 0.0;
        
        for(int k =0; k <ndim; k++) {
          temp += (perm[k]*GradPhi[pres][lp][k][m]*GradPhi[pres][rp][k][m]);
        }
        aux+=( temp * lambdat[m]*JW[m] );
      }
      mx[npos(ntot,mi,mj)]+=aux;
    } // Fim de PE_VI_PP
    
    for(int ls=0;ls<ns;ls++) { // PE_VI_SP
      mj=indice(n0,sat,ls);
      
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      //cout<< "   <- Volume Integrals Ia " << npos(ntot,mi,mj) << "\n";
      for(m=0;m<NGQP;m++) {
        
        double Vaux[ndim];
        for(k=0;k < ndim;k++) {
          //cout<< "   <- sat ls k m "  << sat << "," << ls<< "," << k<< ","<< m <<"\n";
          Vaux[k]=( d_lambdat[m]*Kgrad_pw[k][m]+
                   d_lambdan[m]*Kgrad_pc[k][m]+
                   lambdan[m]*d2_pc[m]*Kgrad_sn[k][m]
                   ) *
          phi_l[m] +
          lambdan[m]*d_pc[m]*perm[k]*GradPhi[sat][ls][k][m];
          //cout<< "Volume Integrals Ia k =" << k << "\n" ;
        }
        
        aaux=0.0;
        
        for(k=0;k<ndim;k++) {
          aaux+=Vaux[k]*GradPhi[pres][rp][k][m];
        }
        
        aux+=(aaux*JW[m]);
        
      }
      
      mx[npos(ntot,mi,mj)]+=aux;
      
    }// Fim de PE_VI_SP
  } // Fim de PE_VI
  
  for(int rs=0;rs<ns;rs++) { // SE_VI
    //   mi=indice(n0,sat,rs);
    // *****************************
    // calculo do vetor
    // *****************************
    
    // primeira integral
    double btemp;
    ptr_stdel[sat]->eval_Phi(rs,phi_r);
    
    aux=0.0;
    for(m=0;m<NGQP;m++) {
      aux+=JW[m]*(sn[m]-sna[m])*phi_r[m];
    }
    aux*=(-porosidade/Dt);
    
    // segunda integral
    aaux=0.0;
    for(m=0;m<NGQP;m++) {
      double temp = 0.0;
      for (k=0; k< ndim; k++) {
        temp += Kgrad_pw[k][m]*GradPhi[sat][rs][k][m];
      }
      aaux+=JW[m]*lambdaw[m]*temp;
      res0[m]=temp;
    }
    aux+=aaux;
    //B[gbnmap[sat][rs]]+=aux; // inicializa valor de B em SE_VI
    btemp=aux;
    // termo de fonte
    
    aaux=qw;
    if(aaux!=0.0) {
      integral(NGQP,JW,phi_r,aux);
      //B[gbnmap[sat][rs]] -= aux*aaux;
      btemp -= aux*aaux;
    }
    //B[gbnmap[sat][rs]] += btemp ;
    mi = indice(n0,sat,rs);
    B[mi] += btemp ;
    // *************************
    // Matriz
    // *************************
    
    for(int lp=0;lp<np;lp++) { // SE_VI_PS
      mj=indice(n0,pres,lp);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        
        double temp =0.0;
        for ( int k=0; k < ndim; k++) {
          temp += (perm[k]*GradPhi[pres][lp][k][m]*GradPhi[sat][rs][k][m]);
        }
        
        aux+=(temp *lambdaw[m]*JW[m]);
      }
      mx[npos(ntot,mi,mj)]+=aux;
    }
    
    for(int ls=0;ls<ns;ls++) { // SE_VI_SS
      mj=indice(n0,sat,ls);
      ptr_stdel[sat]->eval_Phi(ls,phi_l);
      aux=0.0;
      for(m=0;m<NGQP;m++) {
        aux+=JW[m]*(d_lambdaw[m]*res0[m]-porosidade/Dt*phi_r[m])*phi_l[m];
      }
      mx[npos(ntot,mi,mj)]+=aux;
    }
  }
  // ***************************************************************************
  // Colocar os indices no vetor Ti
  // ***************************************************************************
  for(int rp=0;rp<np;rp++) {
    mi=indice(n0,pres,rp);
    indx[mi] = gbnmap[pres][rp];
  }
  for(int rs=0;rs<ns;rs++) {
    mi=indice(n0,sat,rs);
    indx[mi] = gbnmap[sat][rs];
  }
  A->InsertGlobalValues(ntot,indx,mx,Epetra_FECrsMatrix::ROW_MAJOR);
  RHS->SumIntoGlobalValues(ntot,indx,B);
  //  for(int i = 0; i < ntot*ntot; i++) {
  //    printf("mx[%3d] = %12.5e\n",i,mx[i]);
  //  }
  // delete [] indx;indx=nullptr;
  // delete [] mx; mx=nullptr;
  // delete [] B; B=nullptr;
};
// alterado em 21/09/2011
// alterado em 21/10/2011
// alterado em 23/10/2011
// alterado em 13/02/2013

// **************************************************************
// Inicaliza os vetores locais: JV,b0,bs,Grad e traco
// **************************************************************

template<int NumVariaveis>
void PhElem<NumVariaveis>::inicia_vetores()
{
  int nn;//nb,q0,q1;
  int h,i,j,k;
  
  if(vetores_iniciados != 0) { cout<< "vetores locais de PhElem já iniciados\n"; }
  
  else {
    vetores_iniciados = 1;
    //cout<< "PhElem<NumVariaveis>::inicia_vetores()\n";
    for (k=0;k<NumVariaveis;++k){
      nn=ptr_stdel[k]->nn_val();
      u0[k] = new double [nn];
      //  ua[k] = new double [nn];
      usave[k] = new double [nn];
      for(i=0;i<nn;++i){
        u0[k][i]=0.0;
      }
      // nb=ptr_stdel[k]->nb_val();
    }
    // **************************************
    // Calcular o Jacobiano   *             *
    // **************************************
    // q0= ptr_stdel[0]->Q_val(0);
    // q1= ptr_stdel[0]->Q_val(1);
    
    int ndim =  ptr_stdel[0]->ndim_val();
    int qmax = ptr_stdel[0]->qborder_val();
    int nborder = ptr_stdel[0]->nborder_val();
    int NGQP;
    //printf("Calculo do Jacobiano: dimensao=%d\n",NGQP);
    //JV = new double [NGQP]; // opcao 2: um unico vetor
    //ptr_stdel[0]->Jacobian(ptvert,Vert_map,JV);
    compute_JV(0);
    sna = new double [ ptr_stdel[0]->NGQP_val() ];
    pwa = new double [ ptr_stdel[1]->NGQP_val() ];
    
    // ***********************************************************
    // Alocacao dinamica de memoria para a matriz de massa de sn *
    // Mass_sn[i][j]                                             *
    // ***********************************************************
    int ns=ptr_stdel[0]->nn_val();
    Mass_sn = new double * [ns];
    for(i=0;i<ns;i++)
      Mass_sn[i] = new double [ns];
    // ******************
    // Calcular Mass_sn *
    // ******************
    for(i=0;i<ns;i++){
      Mass_sn[i][i]= ptr_stdel[0]->mass(i,i,JV);
      for(j=i+1;j<ns;j++){
        Mass_sn[i][j] = ptr_stdel[0]->mass(i,j,JV);
        Mass_sn[j][i]=Mass_sn[i][j];
      }
    }
    
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // do vetor Jacobiano Jb
    // [lado] [var] [indice] [direcao] [posicao]   *
    // *********************************************
    Jb = new double [nborder*qmax];
    //cout << "INICIOU JB !!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // **************************************
    // Alocacao dinamica de memoria para os *
    // gradientes dos modos                 *
    // [var] [modo] [direcao] [posicao]   *
    // **************************************
    GradPhi = new double *** [NumVariaveis];// variavel i NumVariaveis=2
    for(i=0;i<NumVariaveis;++i){
      nn=ptr_stdel[i]->nn_val();
      NGQP=ptr_stdel[i]->NGQP_val();
      GradPhi[i] = new double ** [nn];// modo j
      for(j=0;j<nn;j++){
        GradPhi[i][j] = new double * [ndim];// direcao k
        for(k=0;k<ndim;k++){
          GradPhi[i][j][k] = new double [NGQP];// posicao
        }
      }
    }
    // **************************************
    // Alocacao dinamica de memoria para os *
    // Laplacianos  dos modos               *
    // [var] [modo] [posicao]             *
    // **************************************
    LaplacianoPhi = new double ** [NumVariaveis];// variavel i  NumVariaveis=2
    for(i=0;i<NumVariaveis;++i){
      nn  =ptr_stdel[i]->nn_val();
      NGQP=ptr_stdel[i]->NGQP_val();
      LaplacianoPhi[i] = new double * [nn];// indice j
      for(j=0;j<nn;j++){
        LaplacianoPhi[i][j] = new double [NGQP];// posicao
      }
    }
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [var] [lado] [modo] [direcao] [posicao]   *
    // *********************************************
    nborder=ptr_stdel[0]->nborder_val();
    TrGradPhi = new double **** [NumVariaveis];// variavel i  NumVariaveis=2
    for(i=0;i<NumVariaveis;++i){
      nn=ptr_stdel[i]->nn_val();
      TrGradPhi[i] = new double *** [nborder];//lado h
      for(h=0;h<nborder;h++){
        TrGradPhi[i][h] = new double ** [nn];// modo j
        for(j=0;j<nn;j++){
          TrGradPhi[i][h][j] = new double * [ndim];// direcao k
          for(k=0;k<ndim;k++){
            TrGradPhi[i][h][j][k] = new double [qmax];// posicao
          }
        }
      }
    }
    
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [var] [lado] [modo] [posicao]               *
    // *********************************************
    TrKgradPhi_n = new double *** [NumVariaveis];//variavel i  NumVariaveis=2
    for(i=0;i<NumVariaveis;++i){
      nn=ptr_stdel[i]->nn_val();
      TrKgradPhi_n[i] = new double ** [nborder];// lado h
      for(h=0;h<nborder;h++){
        TrKgradPhi_n[i][h] = new double * [nn];// modo j
        for(j=0;j<nn;j++){
          TrKgradPhi_n[i][h][j] = new double [qmax];// posicao
        }
      }
    }
    
    // ****************************************************
    // Alocacao dinamica de memoria para os tracos de phi *
    // [var] [lado] [posicao modo]  [posicao q]           *
    // ****************************************************
    TrPhi = new double *** [NumVariaveis];// variavel i  NumVariaveis=2
    for(i=0;i<NumVariaveis;++i){
      TrPhi[i] = new double ** [nborder]; //lado h
      for(h=0;h<nborder;h++){
        int pmax = 1;
        for(int t=0; t < (ndim-1); t++) {
          pmax *= ( (ptr_stdel[i]->P_val(t)) + 1);// pmax = P +1
        }
        TrPhi[i][h] = new double * [pmax];// posicao do modo j
        for(j=0; j < pmax; j++){
          TrPhi[i][h][j] = new double  [qmax];// posicao de q
        }
      }
    }
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [lado][posicao]                             *
    // *********************************************
    Trsn = new double * [nborder];//lado h
    for(h=0;h<nborder;h++){
      Trsn[h] = new double [qmax];// posicao
    }
    
    Trpw = new double * [nborder];//lado h
    for(h=0;h<nborder;h++){
      Trpw[h] = new double [qmax];// posicao
    }
    // *********************************************
    // Alocacao dinamica de memoria para os tracos *
    // [lado] [direcao] [posicao]                  *
    // *********************************************
    TrKgrad_sn = new double ** [nborder];//lado h
    for(h=0;h<nborder;h++){
      TrKgrad_sn[h] = new double * [ndim];// direcao k
      for(k=0;k<ndim;k++){
        TrKgrad_sn[h][k] = new double  [qmax];// posicao
      }
    }
    TrKgrad_pw = new double ** [nborder];//lado h
    for(h=0;h<nborder;h++){
      TrKgrad_pw[h] = new double * [ndim];// direcao k
      for(k=0;k<ndim;k++){
        TrKgrad_pw[h][k] = new double  [qmax];// posicao
      }
    }
    TrKgrad_pc = new double ** [nborder];//lado h
    for(h=0;h<nborder;h++){
      TrKgrad_pc[h] = new double * [ndim];// direcao
      for(k=0;k<ndim;k++){
        TrKgrad_pc[h][k] = new double  [qmax];// posicao
      }
    }
  }
};

// ****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::finaliza_vetores()
{
  if(vetores_iniciados == 1) {
    vetores_iniciados = 0;
    //cout<< "PhElem<NumVariaveis>::finaliza_vetores()\n";
    delete [] JV; JV=nullptr;
    delete [] Jb; Jb=nullptr;
    
    int i,j,k;
    for (k=0;k<2;k++){
      delete [] u0[k]; u0[k]=nullptr;
      //   delete [] ua[k];ua[k]=nullptr;
      delete [] usave[k];usave[k]=nullptr;
    }
    delete [] stgbtrbmapM;stgbtrbmapM=nullptr;
    delete [] stgbtrbmapP;stgbtrbmapP=nullptr;
    
    delete [] sna;sna=nullptr;
    delete [] pwa;pwa=nullptr;
    
    int ns=ptr_stdel[0]->nn_val();
    for(i=0;i<ns;i++){
      delete [] Mass_sn[i]; Mass_sn[i]=nullptr;
    }
    delete [] Mass_sn; Mass_sn=nullptr;
    
    int ndim =ptr_stdel[0]->ndim_val();
    for(i=0;i<2;i++){
      int nn=ptr_stdel[i]->nn_val();
      for(j=0;j<nn;j++){
        for(k=0;k<ndim;k++){
          delete [] GradPhi[i][j][k]; GradPhi[i][j][k]=nullptr;// posicao
        }
        delete [] GradPhi[i][j];GradPhi[i][j]=nullptr;
        delete [] LaplacianoPhi[i][j]; LaplacianoPhi[i][j]=nullptr;
      }
      delete [] GradPhi[i];GradPhi[i]=nullptr;
      delete [] LaplacianoPhi[i]; LaplacianoPhi[i]=nullptr;
    }
    delete [] GradPhi; GradPhi=nullptr;
    delete [] LaplacianoPhi; LaplacianoPhi=nullptr;
    
    // traco de grad phi
    int nborder=ptr_stdel[0]->nborder_val();
    for(int i=0;i<2;i++){
      for(int h=0;h<nborder;h++){
        int nn=ptr_stdel[i]->nn_val();
        for(j=0;j<nn;j++){
          for(k=0;k<ndim;k++){
            delete [] TrGradPhi[i][h][j][k];TrGradPhi[i][h][j][k]=nullptr;// posicao
          }
          delete [] TrGradPhi[i][h][j];TrGradPhi[i][h][j]=nullptr;
        }
        delete [] TrGradPhi[i][h];TrGradPhi[i][h]=nullptr;
      }
      delete [] TrGradPhi[i]; TrGradPhi[i]=nullptr;
    }
    delete [] TrGradPhi;TrGradPhi=nullptr;
    // traco de grad phi dot n
    for(int i=0;i<2;i++){
      for(int h=0;h<nborder;h++){
        int nn=ptr_stdel[i]->nn_val();
        for(j=0;j<nn;j++){
          delete [] TrKgradPhi_n[i][h][j];TrKgradPhi_n[i][h][j]=nullptr;
        }
        delete [] TrKgradPhi_n[i][h];TrKgradPhi_n[i][h]=nullptr;
      }
      delete [] TrKgradPhi_n[i];TrKgradPhi_n[i]=nullptr;
    }
    delete [] TrKgradPhi_n;TrKgradPhi_n=nullptr;
    
    // traco de Phi
    for(int i=0;i < 2;i++){
      int p0=ptr_stdel[i]->P_val(0);
      for(int h=0;h < nborder;h++){
        for(j=0;j<p0;j++){
          delete [] TrPhi[i][h][j]; TrPhi[i][h][j]=nullptr;
        }
        delete [] TrPhi[i][h];TrPhi[i][h]=nullptr;
      }
      delete [] TrPhi[i];TrPhi[i]=nullptr;
    }
    delete [] TrPhi;TrPhi=nullptr;
    // *********************************************
    // Liberacao de memoria para os tracos         *
    // [lado][posicao]                             *
    // *********************************************
    int h;
    for(h=0;h<nborder;h++){
      delete [] Trsn[h];Trsn[h]=nullptr;
      delete [] Trpw[h];Trpw[h]=nullptr;
    }
    delete [] Trsn;Trsn=nullptr;
    delete [] Trpw;Trpw=nullptr;
    for(h=0;h<nborder;h++){
      for(k=0;k<ndim;k++){
        delete [] TrKgrad_sn[h][k];TrKgrad_sn[h][k]=nullptr;
        delete [] TrKgrad_pw[h][k];TrKgrad_pw[h][k]=nullptr;
        delete [] TrKgrad_pc[h][k];TrKgrad_pc[h][k]=nullptr;
      }
      delete [] TrKgrad_sn[h]; TrKgrad_sn[h]=nullptr;
      delete [] TrKgrad_pw[h]; TrKgrad_pw[h]=nullptr;
      delete [] TrKgrad_pc[h]; TrKgrad_pc[h]=nullptr;
    }
    delete [] TrKgrad_sn;TrKgrad_sn=nullptr;
    delete [] TrKgrad_pw;TrKgrad_pw=nullptr;
    delete [] TrKgrad_pc;TrKgrad_pc=nullptr;
  }
};

// ***************************************************************************
// Salva as rotinas abaixo para servir de template para casos de
// inicializacao de vetores mais simples que o DG_PhElem
/*
 void PhElem::inicia_vetores()
 {
 for (int var=0; var < NumVariaveis; ++var)
 {
 int Nn = pt_stdel[var]->nn_val();
 b0[var] = new double [Nn];
 u0[var] = new double [Nn];
 //u1 = new double [Nn];
 //u2 = new double [Nn];
 for(int i=0;i<Nn;i++){
 b0[var][i]=0.0;
 u0[var][i]=0.0;
 // u1[i]=0.0;
 //u2[i]=0.0;
 }
 int Nb=pt_stdel[var]->nb_val();
 bs[var] = new double [Nb];
 for(int i=0;i<Nb;i++)bs[var][i]=0.0;
 }
 };
 // *************************************************************
 void PhElem::finaliza_vetores()
 {
 for (int var=0; var < NumVariaveis; ++var)
 {
 delete [] b0[var];  b0[var] = nullptr;
 delete [] u0[var];  u0[var] = nullptr;
 delete [] bs[var];  bs[var] = nullptr;
 }
 delete [] b0; b0 = nullptr;
 delete [] u0; u0 = nullptr;
 delete [] bs; bs = nullptr;
 };
 
 // *************************************************************
 void PhElem::zera_vetores()
 {
 for (int var=0;var < NumVariaveis; ++var)
 {
 int Nn=pt_stdel[var]->nn_val();
 for(int i=0;i<Nn;i++){
 b0[var][i]=0.0;
 u0[var][i]=0.0;
 //u1[i]=0.0;
 //u2[i]=0.0;
 }
 int Nb=pt_stdel[var]->nb_val();
 for(int i=0;i<Nb;i++)bs[var][i]=0.0;
 }
 };
 void PhElem::zera_bs()
 {
 for (int var=0;var < NumVariaveis; ++var)
 {
 int Nb=pt_stdel[var]->nb_val();
 for(int i=0;i<Nb;i++)bs[var][i]=0.0;
 }
 };
 */

// *****************************************************************************
//#include "GBNMAP.cpp"

// *****************************************************************************
// Mapeamento dos vertices
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_vertices(std::vector< std::vector<int>> gbn)//int ** gbn )
{
  for(int i=0;i<NumVariaveis;++i){
    for(int j=0;j<numv;++j){
      gbnmap[i][j]=gbn[i][Vert_map[j]];
      sgn[i][j]=1;
    }
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_vertices(int gbn[],const int & ivar)
{
  for(int i=0;i<numv;++i){
    gbnmap[ivar][i]=gbn[Vert_map[i]];
    sgn[ivar][i]=1;
  }
};

// *****************************************************************************
// Mapeamento das arestas
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_aresta(const int & aresta,const int & sinal,
                                               const int & ivar,const int & inicio,
                                               int & inc)
{
  int ip;
  int ini = ptr_stdel[ivar]->show_emapi(aresta)+1;
  int fim = ptr_stdel[ivar]->show_emapi(aresta+1)-1;
  int aux=1;
  //cout << "(inicio,fim) = (" << ini << ","<< fim << ")   inc_inicial " << inc << endl;
  for(int i = ini; i < fim; ++i){
    ip=ptr_stdel[ivar]->show_emapv(i);
    gbnmap[ivar][ip]=inicio+inc;
    sgn[ivar][ip]=aux;
    aux*=sinal; // muda o sinal de aux se sinal=-1
    ++inc;
  }
  //cout << "inc_final "<< inc << endl;
};

// *****************************************************************************
// Mapeamento dos modos interiores
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_interior(const int & ivar, int & count)
{
  int ini=ptr_stdel[ivar]->nb_val();
  int fim=ptr_stdel[ivar]->nn_val();
  for(int i=ini;i<fim;++i){
    gbnmap[ivar][i]=count++;
    sgn[ivar][i]=1;
  }
};

// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_faces(const int face_vec[],const int & ivar)
{
  int i[3];
  cout << "PhElem<NumVariaveis>::gbnmap_faces\n";
  // Loop sobre as faces
  for(int f=0;f<numf;++f){
    
    int dir[3];
    int v2=0;
    int sinal[2];
    
    // achar os vertices da face
    int nvf = ptr_stdel[ivar]->show_nvf(f);
    int Av[nvf], ind[nvf];
    for(int i = 0; i < nvf; ++i){
      ind[i]=ptr_stdel[ivar]->face_lvert(f,i);
      Av[i]=Vert_map[ind[i]];
    }
    
    if(nvf == 4) {
      quad_ordem(Av,ind,dir,sinal,v2);
      
      if(v2 == 0) i[dir[2]] = 0;
      else i[dir[2]] = ptr_stdel[ivar]->P_val(dir[2]);
      int P0 = ptr_stdel[ivar]->P_val(dir[0]);
      int P1 = ptr_stdel[ivar]->P_val(dir[1]);
      
      int inicio = face_vec[face_map[f]];
      int inc = 0;
      int aux0 = 1;
      
      for(int i0 = 1; i0 < P0; ++i0) {
        i[dir[0]] = i0;
        int aux1 = 1;
        for(int i1 = 1; i1< P1; ++i1) {
          i[dir[1]] = i1;
          int n = ptr_stdel[ivar]->show_ind_mode(i[0],i[1],i[2]);
          gbnmap[ivar][n]=inicio + inc;
          sgn[ivar][n]=aux1*aux0;
          inc++;
          aux1 *= sinal[1];
        }
        aux0 *= sinal[0];
      }
    }
    else if ( nvf == 3) { // faces triangulares
      dir[0]=ptr_stdel[ivar]->show_fd0(f);
      dir[1]=ptr_stdel[ivar]->show_fd1(f);
      dir[2]=3-dir[0]-dir[1];
      v2=ptr_stdel[ivar]->show_fv2(f);
      
      if(v2 == 0) i[dir[2]] = 0;
      else i[dir[2]] = ptr_stdel[ivar]->P_val(dir[2]);
      int P0 = ptr_stdel[ivar]->P_val(dir[0]);
      int P1 = ptr_stdel[ivar]->P_val(dir[1]);
      
      int inicio = face_vec[face_map[f]];
      int inc = 0;
      //int aux0 = 1;
      
      for(int i0 = 1; i0 < P0; ++i0) {
        i[dir[0]] = i0;
        //int aux1 = 1;
        for(int i1 = 1; i1< P1 - i0; ++i1) {
          i[dir[1]] = i1;
          int n = ptr_stdel[ivar]->show_ind_mode(i[0],i[1],i[2]);
          gbnmap[ivar][n]=inicio + inc;
          sgn[ivar][n]=1;//aux1*aux0;
          inc++;
          //aux1 *= sinal[1];
        }
        //aux0 *= sinal[0];
      }
    }
  }
};
// *****************************************************************************
template<int NumVariaveis>
void PhElem<NumVariaveis>::gbnmap_arestas(const int aresta_vec[],const int & ivar)
{
  int av[2];
  // Loop sobre as arestas
  for(int a=0;a<nume;++a){
    av[0]=ptr_stdel[ivar]->aresta_lvert(a,0);
    av[1]=ptr_stdel[ivar]->aresta_lvert(a,1);
    int sinal = 1;
    if(Vert_map[av[0]] > Vert_map[av[1]]) sinal = -1;
    
    int ini = ptr_stdel[ivar]->show_emapi(a)+1;
    int fim = ptr_stdel[ivar]->show_emapi(a+1)-1;
    int inicio = aresta_vec[aresta_map[a]];
    int aux=1;
    int inc = 0;
    for(int i = ini; i < fim; ++i){
      int ip=ptr_stdel[ivar]->show_emapv(i);
      gbnmap[ivar][ip]=inicio+inc;
      sgn[ivar][ip]=aux;
      aux*=sinal; // muda o sinal de aux se sinal=-1
      inc++;
    }
  }
};

#endif /* PhElem_hpp */
