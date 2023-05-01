import torch
import numpy as np
from torch import nn

grid=100

xgr=torch.asarray([[(i-1)/(grid-1)]for i in range(1,grid+1)]).requires_grad_(True)
ugr=torch.asarray([[2*(i-1)/(grid-1)-1]for i in range(1,grid+1)]).requires_grad_(True)


class Neury(nn.Module):
    def __init__(self):
        super().__init__()  #call the superclass nn.Module with which arguments??
        #self.flatten = nn.Flatten() #TAKES THE TENSOR ARRAY INTO A LIST?
        self.radial = nn.Sequential(
            nn.Linear(1, 200),  #the input in 28*28
            nn.Tanh(),
            nn.Linear(200, 2),
            nn.Tanh(),
        )
        self.angular = nn.Sequential(
            nn.Linear(1, 200),  #the input in 28*28
            nn.Tanh(),
            nn.Linear(200, 2),
            nn.Tanh(),
        ) 

    def hard_enforcement_rad(self,x,p):
        return (torch.mul(torch.exp(x-1),p))+1 #torch.mul may be just multiplication?
    
    def hard_enforcement_ang(self,x,p):
        return (torch.mul(torch.exp(x+1),p))+1 


    #the pytorch class Module calls a method "forward" automatically.
    def forward(self, x,y):
        outx = self.radial(x)
        frad = self.hard_enforcement_rad(x,outx) #maybe this class method should be defined as function outside the class
        outy = self.angular(y)
        gang = self.hard_enforcement_ang(y,outy)
        return [frad,gang]
PINN = Neury()

#########################################
#coding the ODEs. Check the coefficients!
ell=2
ss=0
rp=1
m=0
freq=torch.as_tensor(0.7 - 0.1j).requires_grad_(True)
A=ell*(ell+1)- ss*(ss+1)

def F1(x,w):
    return (A*(x-1) + (1 + ss)*x + rp*(1 + ss)*(-2 + x)*x - 4*rp*w*w*(-1 + x**2))+2*w*1j*(-1+rp-ss+2*rp*(ss+1)*x - rp*(2+ss)*x**2)
def F2(x,w):
    return (x-1)*(2*rp*x*x+(ss+1)*(x*x-2*x))+1j*(x-1)*(-4*rp*w+2*w)
def F3(x,w):
    return (x**2)*((x-1)**2) + 1j*0
def G1(u,w):
    return 4*(A*(u**2-1)+m**2+2*m*ss*u+ss*((ss+1)*u**2-1))-2*(u**2-1)*abs(m-ss)-2*(u**2)*abs(m-ss)*(abs(m+ss)+1)-((m-ss)**2)*(u-1)**2-((u+1)**2)*(m+ss)**2
def G2(u,w):
    return -4*(u**2-1)*((u-1)*abs(m-ss)+(u+1)*abs(m+s)+2*u)
def G3(u,w):
    return -4*(u**2-1)**2

def dfx(x,f):
    return torch.autograd.grad([f], x, grad_outputs=torch.ones(x.squeeze().shape), create_graph=True)[0]
##############################################




optimizer= torch.optim.Adam(PINN.parameters(),lr=0.005,weight_decay=0.999)
optimizer.add_param_group({'params':freq})  #optimizer.param_groups.append({'params': extra_params })



PINN.train()
for epoch in range(2000):

    xvar=xgr.squeeze()
    uvar=ugr.squeeze()

    fr,fi=PINN(xgr,ugr)[0].t()
    gr,gi=PINN(xgr,ugr)[1].t()

    frp=dfx(xgr,fr).squeeze()
    frpp=dfx(xgr,frp).squeeze()
    fip=dfx(xgr,fi).squeeze()
    fipp=dfx(xgr,fip).squeeze()
    grp=dfx(ugr,gr).squeeze()
    grpp=dfx(ugr,grp).squeeze()
    gip=dfx(ugr,gi).squeeze()
    gipp=dfx(ugr,gip).squeeze()

    fc=torch.view_as_complex(torch.stack((fr,fi)).t().contiguous())
    fcp=torch.view_as_complex(torch.stack((frp,fip)).t().contiguous())
    fcpp=torch.view_as_complex(torch.stack((frpp,fipp)).t().contiguous())
    gc=torch.view_as_complex(torch.stack((gr,gi)).t().contiguous())
    gcp=torch.view_as_complex(torch.stack((grp,gip)).t().contiguous())
    gcpp=torch.view_as_complex(torch.stack((grpp,gipp)).t().contiguous())

    WW=10
    losss=WW*sum(torch.abs((F1(xvar,freq)*fc+F2(xvar,freq)*fcp+F3(xvar,freq)*fcpp)))+sum(torch.abs((G1(uvar,freq)*gc+F2(uvar,freq)*gcp+F3(uvar,freq)*gcpp)))

    optimizer.zero_grad()
    losss.backward()
    optimizer.step()
    print(losss,freq)








