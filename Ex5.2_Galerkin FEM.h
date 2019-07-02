load "UMFPACK64"
int NN=16; 
int NN2=1024;
real dt=1.0/NN, t=0.0, w=3.0*pi/2.0, beta=pi/w; 
complex I=1.0i; 

//Define the boundary of the computational domain
border aa(t=0,0.5){x=0; y=-t; label=1;}; 
border ab(t=0,0.5){x=t; y=-0.5; label=2;}; 
border ac(t=0,0.5){x=0.5; y=t-0.5; label=1;}; 
border ad(t=0,0.5){x=0.5-t; y=0; label=2;}; 
border ba(t=0,1.5){x=-0.5; y=0.5-t; label=3;}; 
border bb(t=0,1.5){x=t-0.5; y=-1.0; label=4;}; 
border bc(t=0,1.5){x=1.0; y=t-1.0; label=3;}; 
border bd(t=0,1.5){x=1.0-t; y=0.5; label=4;}; 

//Triangulation, with NN nodes per unit length on the boundary 
mesh Th = buildmesh (ba(3*NN/2)+bb(3*NN/2)+bc(3*NN/2)+bd(3*NN/2)
+aa(-NN/2)+ab(-NN/2)+ac(-NN/2)+ad(-NN/2));   
//plot(Th,ps="Th.eps",wait=true); 
mesh Th2 = buildmesh (ba(3*NN2/2)+bb(3*NN2/2)+bc(3*NN2/2)+bd(3*NN2/2)
+aa(-NN2/2)+ab(-NN2/2)+ac(-NN2/2)+ad(-NN2/2)); 

//Define the P2-P1 Taylor-Hood finite element space
fespace Sh(Th,P2); 
fespace Vh(Th,P1);  
Sh u1,u2,H1,H2,v1,v2,w1,w2,u10,u20,H10,H20; 
Vh p,q; 
Sh varphiR;

//Compute a surrogate for the exact harmonic function varphi by using a much smaller mesh
fespace Sh2(Th2,P2);
Sh2 varphi2,phi2;  
problem VPHI2(varphi2,phi2,solver=sparsesolver) //solver=GMRES
=int2d(Th2)(dx(varphi2)*dx(phi2)+dy(varphi2)*dy(phi2))      
+on(1,2,varphi2=1.0)+on(3,4,varphi2=0.0);  
VPHI2;
varphiR=varphi2;

//Initial condition
u10=0.0; 
u20=0.0; 
H10=dy(varphiR); 
H20=-dx(varphiR);  

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
func curlH=-D2u1;
func f1=u1t+(u1e*u1x+u2e*u1y)-D2u1+px+H2e*curlH;
func f2=u2t+(u1e*u2x+u2e*u2y)-D2u2+py-H1e*curlH;
func J=u1t+curlH-u1e*(H2e-H1e);
func g=u1x+u2y;
 
//Weak formulation of the standard Galerkin FEM for MHD
problem BackEuler([H1,H2,u1,u2,p],[w1,w2,v1,v2,q],solver=sparsesolver) 
=int2d(Th)((H1*w1+H2*w2)/dt+(u1*v1+u2*v2)/dt
+(dx(H2)-dy(H1))*(dx(w2)-dy(w1))
-(u1*H20-u2*H10)*(dx(w2)-dy(w1))
+(dx(u1)*dx(v1)+dy(u1)*dy(v1) + dx(u2)*dx(v2)+ dy(u2)*dy(v2))
+0.5*(u10*dx(u1)+u20*dy(u1))*v1+0.5*(u10*dx(u2)+u20*dy(u2))*v2
-0.5*(u10*dx(v1)+u20*dy(v1))*u1-0.5*(u10*dx(v2)+u20*dy(v2))*u2
-p*dx(v1)-p*dy(v2) 
+H20*(dx(H2)-dy(H1))*v1-H10*(dx(H2)-dy(H1))*v2
+q*(dx(u1)+dy(u2))  
+1.0e-8*p*q)       
-int2d(Th)((H10*w1+H20*w2)/dt+(J+u1e*(dx(varphiR)+dy(varphiR)))*(dx(w2)-dy(w1))
+(f1-dx(varphiR)*curlH)*v1+(f2-dy(varphiR)*curlH)*v2+(u10*v1+u20*v2)/dt+g*q)  
+on(1,3,u1=0.0,u2=0.0,H1=0.0)+on(2,4,u1=0.0,u2=0.0,H2=0.0);  

//Time stepping
for (int j=0; j<1.0/dt; j=j+1) 
{ 
t=t+dt;
BackEuler;
H10=H1;
H20=H2;
u10=u1;
u20=u2;
cout<<"t="<<t<<endl;
} 

//Display the error of the numerical solution
cout<<"L^2 error of H="<<sqrt(int2d(Th)(abs(H1e+dy(varphiR)-H1)^2
+abs(H2e-dx(varphiR)-H2)^2))<<endl;
cout<<"L^2 error of u="<<sqrt(int2d(Th)((u1-u1e)^2+(u2-u2e)^2)) <<endl; 
