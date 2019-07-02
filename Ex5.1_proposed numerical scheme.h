int NN=16; 
real dt=1.0/NN, t=0.0, w=3.0*pi/2.0, beta=pi/w; 
complex I=1.0i; 

//Define the boundary of the computational domain
border ba(t=0,0.5){x=t; y=0.0; label=1;}; 
border bb(t=0,0.5){x=0.5; y=t; label=1;}; 
border bc(t=0,1.0){x=0.5-t; y=0.5; label=1;}; 
border bd(t=0,1.0){x=-0.5; y=0.5-t; label=1;}; 
border be(t=0,0.5){x=t-0.5; y=-0.5; label=1;}; 
border bf(t=0,0.5){x=0.0; y=t-0.5; label=1;}; 

//Triangulation, with NN nodes per unit length on the boundary
mesh Th = buildmesh (ba(NN/2)+bb(NN/2)+bc(NN)+bd(NN)+be(NN/2)+bf(NN/2) ); 
//plot(Th,ps="Th.eps",wait=true); 

//Define the P2-P1 Taylor-Hood finite element space
fespace Sh(Th,P2); 
fespace Vh(Th,P1);  
Sh u1,u2,A,B,v1,v2,a,wh,u10,u20,A0; 
Vh p,q; 

//Initial condition
u10=0.0; 
u20=0.0; 
A0=0.0; 

//Define the exact solution and right-hand sides
func real HS(real s) { return s>0; }; 
func r=sqrt(x^2+y^2); 
func theta=HS(arg(x+y*I))*arg(x+y*I)+HS(-arg(x+y*I)-1e-10)*(arg(x+y*I)+2.0*pi); 
func ss=(r-0.1)/(0.4-0.1); 
func Phi=1e-1*(1.0+20.0*ss^7-70.0*ss^6+84.0*ss^5-35.0*ss^4)*HS(0.4 -r)*HS(r-0.1)+1e-1*HS(0.1-r); 
func Phi1=1e-1*(140.0*ss^6-420.0*ss^5+420.0*ss^4-140.0*ss^3)*10.0/3.0*HS(0.4 -r)*HS(r-0.1); 
func Phi2=1e-1*(840.0*ss^5-2100.0*ss^4+1680.0*ss^3-420.0*ss^2)*100.0/9.0*HS(0.4 -r)*HS(r-0.1);  
func Phi3=1e-1*(4200.0*ss^4-8400.0*ss^3+5040.0*ss^2-840.0*ss)*1000.0/27.0*HS(0.4 -r)*HS(r-0.1);  
func PPh=(Phi1*(2.0*beta+1.0)*r^(beta-1)+Phi2*r^beta); 
func u1e=t^2*Phi*r^beta*sin(beta*theta); 
func u2e=t^2*Phi*r^beta*sin(beta*theta); 
func u1t=2.0*t*Phi*r^beta*sin(beta*theta); 
func u2t=2.0*t*Phi*r^beta*sin(beta*theta); 
func u1x=t^2*(Phi1*r^beta*sin(beta*theta)*cos(theta)+Phi*beta*r^(beta-1.0)*sin((beta-1.0)*theta)); 
func u1y=t^2*(Phi1*r^beta*sin(beta*theta)*sin(theta)+Phi*beta*r^(beta-1.0)*cos((beta-1.0)*theta)); 
func u2x=t^2*(Phi1*r^beta*sin(beta*theta)*cos(theta)+Phi*beta*r^(beta-1.0)*sin((beta-1.0)*theta)); 
func u2y=t^2*(Phi1*r^beta*sin(beta*theta)*sin(theta)+Phi*beta*r^(beta-1.0)*cos((beta-1.0)*theta)); 
func D2u1=t^2*PPh*sin(beta*theta); 
func D2u2=t^2*PPh*sin(beta*theta); 
func pe=0.0;
func px=0.0;
func py=0.0;
func H1e=u1y; 
func H2e=-u1x;
func Ae=u1;
func curlH=-D2u1;
func f1=u1t+(u1e*u1x+u2e*u1y)-D2u1+px+H2e*curlH;
func f2=u2t+(u1e*u2x+u2e*u2y)-D2u2+py-H1e*curlH;
func J=u1t+curlH-u1e*(H2e-H1e);
func g=u1x+u2y;
 
//Weak formulation of the proposed numerical scheme for MHD
problem BackEuler([A,u1,u2,B,p],[a,v1,v2,wh,q],solver=UMFPACK) 
=int2d(Th)(A*a/dt+(u1*v1+u2*v2)/dt
+(dx(A)*dx(a)+dy(A)*dy(a))
+(u1*dx(A0)+u2*dy(A0))*a
+(dx(u1)*dx(v1)+dy(u1)*dy(v1) + dx(u2)*dx(v2)+ dy(u2)*dy(v2))
+0.5*(u10*dx(u1)+u20*dy(u1))*v1+0.5*(u10*dx(u2)+u20*dy(u2))*v2
-0.5*(u10*dx(v1)+u20*dy(v1))*u1-0.5*(u10*dx(v2)+u20*dy(v2))*u2
-p*dx(v1)-p*dy(v2) 
+dx(A0)*B*v1+dy(A0)*B*v2
+B*wh+dx(A)*dx(wh)+dy(A)*dy(wh) 
+q*(dx(u1)+dy(u2))  
+1.0e-8*p*q)       
-int2d(Th)((A0/dt+J)*a+f1*v1+f2*v2+(u10*v1+u20*v2)/dt+g*q)  
+on(1,A=0.0,u1=0.0,u2=0.0,B=0.0);  

//Time stepping
for (int j=0; j<1.0/dt; j=j+1) 
{ 
t=t+dt;
BackEuler;
A0=A;
u10=u1;
u20=u2;
cout<<"t="<<t<<endl;
} 

//Display the error of the numerical solution
cout<<"L^2 error of A="<<sqrt(int2d(Th)(abs(Ae-A)^2))<<endl;
cout<<"L^2 error of H="<<sqrt(int2d(Th)(abs(H1e-dy(A))^2+abs(H2e+dx(A))^2))<<endl;
cout<<"L^2 error of u="<<sqrt(int2d(Th)((u1-u1e)^2+(u2-u2e)^2)) <<endl; 
