//dynamic array structs

#include <iostream>
using namespace std;

struct flexArrayStruct{
  int num;
  int data[];
  int tryd;
};

/* ... */
int main() {
   size_t array_size = 10;
   size_t i;
   
   /* Initialize array_size */
   
   /* Space is allocated for the struct */
   struct flexArrayStruct *structP = (struct flexArrayStruct *)
      malloc(sizeof(struct flexArrayStruct) + sizeof(int) * array_size);
   if (structP == NULL) {
   /* Handle malloc failure */
      cout<<"malloc failure"<<endl;
   }
   
   structP->num = 0;
   structP->tryd = 0;
   
   /* Access data[] as if it had been allocated
   * as data[array_size]
   */
   for (i = 0; i < array_size; i++) {
   structP->data[i] = i+1;
   }
   
   structP->tryd = 999;
   
   cout<<structP->num<<endl;
   for (i = 0; i < array_size; i++) {
      cout<<structP->data[i]<<", ";
   }
   cout<<endl;
   cout<<structP->tryd<<endl;
}
