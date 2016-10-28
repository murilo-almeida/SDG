#include "spectral.h"
// ****************************************************************************
/*
 //template<int T_VAR>
void  PhElem::FindVarsInRow(int * NumNz, vector<int>  MapRow[])
{
 
  const int qmax=ptr_stdel[sat]->qborder_val();
  const int NGQP=ptr_stdel[sat]->NGQP_val();
  
  const int np=ptr_stdel[pres]->nn_val();
  const int ns=ptr_stdel[sat ]->nn_val();
  int row, col;
  for(int rs=0;rs<ns;rs++) { // SE_VI
    row=gbnmap[sat][rs];
    // *************************
    // Matriz
    // *************************
    for(int ls=0;ls<ns;ls++) { // SE_VI_SS
      col=gbnmap[sat][ls];
      //MapRow[row][NumNz[row]]=col;
 MapRow[row].pushback(col);
      NumNz[row]++;
    }
    for(int lp=0;lp<np;lp++) { // SE_VI_PS
      col=gbnmap[pres][lp];
      //MapRow[row][NumNz[row]]=col;
  MapRow[row].pushback(col);
      NumNz[row]++;
    }
  }
  for(int rp=0;rp<np;rp++) { // PE_VI
    row=gbnmap[pres][rp];
    // *****************************
    // Matriz
    // *****************************
    for(int ls=0;ls<ns;ls++) { // PE_VI_SP
      col=gbnmap[sat][ls];
      //MapRow[row][NumNz[row]]=col;
  MapRow[row].pushback(col);
 NumNz[row]++;
    }// Fim de PE_VI_SP
    for(int lp=0;lp<np;lp++) { // PE_VI_PP
      col=gbnmap[pres][lp];
      //MapRow[row][NumNz[row]]=col;
  MapRow[row].pushback(col);
 NumNz[row]++;
    } // Fim de PE_VI_PP
  } // Fim de PE_VI
}
*/

 /*
  //template<int T_VAR>
void PhElem::Row(int * NumNz, vector<int> * MapRow)
{
  const int NN = 2;
  int indr,indl;
  int nn[NN];
  for(int i=0;i<NN;i++)
    nn[i]=ptr_stdel[i]->nn_val();
  int row, col;
  for(indr=0;indr<NN;indr++) {
    
    for(int r=0;r<nn[indr];r++) { // SE_VI
      row=gbnmap[indr][r];
      // *************************
      // Matriz
      // *************************
      for(indl=0;indl<NN;indl++) {
	
	for(int l=0;l<nn[indl];l++) { // SE_VI_SS
	  col=gbnmap[indl][l];
	  //MapRow[row][NumNz[row]]=col;
	 MapRow[row].pushback(col);
  NumNz[row]++;
	}
      
      }
    
    }
  }
};
 */
// *************************************************************************************************
/*
void DG_Prob::FindVarsInRowEI(const EDGE border, int * NumNz, vector<int> * MapRow)
{
  
  int pos,iEr,iEl;
  int row,col;
  //static int index=0;
  //const int tipo=border.tipo;
  // elementos vizinhos ao border
  PhElem e[2];
  e[0]=el[border.elemento[0]];
  e[1]=el[border.elemento[1]];
  
  // Variaveis para dimensionar os arrays abaixo
 // const int qmax=e[0].show_ptr_stdel(0)->qborder_val();// Pontos de gauss nas arestas
  //const int nsat=e[0].show_ptr_stdel(sat)->nn_val();// Num de modos para saturacao
  //const int npres=e[0].show_ptr_stdel(pres)->nn_val();// Num de modos para pressao
  
  int ns[2],np[2];
  np[0]=e[0].show_ptr_stdel(pres)->nn_val();
  ns[0]=e[0].show_ptr_stdel(sat)->nn_val();
  np[1]=e[1].show_ptr_stdel(pres)->nn_val();
  ns[1]=e[1].show_ptr_stdel(sat)->nn_val();
  int a[2];
  a[0]=border.num_local[0];
  a[1]=border.num_local[1];
  
  
  // **************************************************************************  
  // ********
  // Matriz *
  // ********
  
  for(iEr=0;iEr<2;iEr++) {// iEr                                             // PE_IE_P   
    iEl = (iEr==0) ? 1 : 0; 		    
    // PE_IE_P
    // Pressure								    // PE_IE_P
    // PE_IE_P
    for(int rp=0;rp<np[iEr];rp++) {// PE_IE  rp				    // PE_IE_P
      
      int rp_on_border=e[iEr].show_ptr_stdel(pres)->is_on_border(rp,a[iEr],pos);	    // PE_IE_P
      if(rp_on_border) 
	row=e[iEr].map(pres,rp);
      
      for(int lp=0;lp<np[iEl];lp++) { // PE_IE_P lp			    // PE_IE_P
	
	int lp_on_border=e[iEl].show_ptr_stdel(pres)->is_on_border(lp,a[iEl],pos);    // PE_IE_P
	if(rp_on_border && lp_on_border) {
	  //MapRow[row][NumNz[row]]
 col = e[iEl].map(pres,lp);
 MapRow[row].pushback(col);
 NumNz[row]++;
	}								    // PE_IE_P
      }                                                                   // PE_IE_P 
      
      // Saturation
      
      for(int ls=0;ls<ns[iEl];ls++) { // PE_IE_S ls   		     
	int ls_on_border=e[iEl].show_ptr_stdel(sat)->is_on_border(ls,a[iEl],pos);// PE_IE_S	
	if(rp_on_border && ls_on_border) {
	  //MapRow[row][NumNz[row]]
 col=e[iEl].map(sat,ls);
  MapRow[row].pushback(col);
	  NumNz[row]++;							    
	}							       // PE_IE_S
      }							       // PE_IE_S
    }								       // PE_IE_S
 								       // PE_IE_S
    // PE_IE_S      
    
    // **************************************************************************
    // Saturation Equation
    // **************************************************************************
    // SE_IE_P 
    // SE_IE_P
    for(int rs=0;rs<ns[iEr];rs++){// SE_IE  rs				   // SE_IE_P
      // SE_IE_P
      int rs_on_border=e[iEr].show_ptr_stdel(sat)->is_on_border(rs,a[iEr],pos);	   // SE_IE_P
      if(rs_on_border)
	row=e[iEr].map(sat,rs);
      // SE_IE_P
      for(int lp=0;lp<np[iEl];lp++){ 					 
	int lp_on_border=e[iEl].show_ptr_stdel(pres)->is_on_border(lp,a[iEl],pos);   // SE_IE_P
	if(rs_on_border && lp_on_border) {
	  //MapRow[row][NumNz[row]]
 col=e[iEl].map(pres,lp);
 MapRow[row].pushback(col);
 NumNz[row]++;
	}								
      }                                                                 
      
      for(int ls=0;ls<ns[iEl];ls++) {              						    
	int ls_on_border=e[iEl].show_ptr_stdel(sat)->is_on_border(ls,a[iEl],pos);
	if(rs_on_border && ls_on_border) {
	  //MapRow[row][NumNz[row]]
 col=e[iEl].map(sat,ls);
 MapRow[row].pushback(col);
 NumNz[row]++;
	}							   
      }	
    } 
  } 
};
*/
