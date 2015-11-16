#include "makeWorkspace.C"

void fit2015() {
   RooWorkspace myws;
   makeWorkspace(myws, 1, -2.4, 2.4, 0, 100, 0, 2.4, 0, 50, 0, 200);

   buildModelUpsi2015(myws, 2, 3);

}
