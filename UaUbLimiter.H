if(InitFreeze.value()==1)
{
// in this way the bed is immobile in streamwise
    volScalarField temp = Ua & g;
    volScalarField temp2 = Ub & g;
// extract vertical component of velocity component,we get ( 0 ,0,w)
    volVectorField Uatemp = temp*g/magSqr(g);
    volVectorField Ubtemp = temp2*g/magSqr(g);
// extract spanwise component of velocity, we get ( 0,v,0)
    volVectorField Uatemp2 = (Ua & vector(0,1,0))*vector(0,1,0);

    scalar indv1 = (vector(0,1,0) & g.value());
    scalar indv2 = (vector(0,0,1) & g.value());

    int ind = 0;
    if(indv1 !=0) ind = 1;
    if(indv2 !=0) ind = 2;

    forAll (alpha, celll)
    {
//       if(alpha[celll]>=alphaMinFriction.value())
//       {
// in this way, the fluid veocity at the bottom would not be affected by the wall function.
           Ua[celll] = Uatemp[celll];
           Ub[celll] = Ubtemp[celll];
//       }
    }
    Ua.correctBoundaryConditions();
    Ub.correctBoundaryConditions();
}
if(InitFreeze.value()==2)
{
// in this way the bed is immobile in streamwise
   volScalarField temp = Ua & g;
// extract vertical component of velocity component,we get ( 0 ,0,w)
   volVectorField Uatemp = temp*g/magSqr(g);
// extract spanwise component of velocity, we get ( 0,v,0)
   volVectorField Uatemp2 = (Ua & vector(0,1,0))*vector(0,1,0);
   scalar indv1 = (vector(0,1,0) & g.value());
    scalar indv2 = (vector(0,0,1) & g.value());

    int ind = 0;
    if(indv1 !=0) ind = 1;
    if(indv2 !=0) ind = 2;
    
    forAll (alpha, celll)
    {
 //      if(alpha[celll]>=alphaMinFriction.value())
 //      {
          Ua[celll] = Uatemp[celll]+Uatemp2[celll];
 //      }
    }
    Ua.correctBoundaryConditions();
}
