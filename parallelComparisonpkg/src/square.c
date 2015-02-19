//#include <stdio.h>

//int main(){
//  return 0;
//}

void square(int *N , double *x){
  /*Takes an array of size N and squares all entries.*/
  int n = *N;
  //size_t i = (size_t) n;
  int i = (int) n;
  for(i=0;i<n;++i){
    x[i] = x[i] * x[i];
  }
}
