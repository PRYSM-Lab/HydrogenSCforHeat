****************************************************************
*Hydrogen Infrastructure
*Electrolysis from renewables or grid
*HA2
*HA1
****************************************************************
$set GAP 0.05
$set outputname monolithic_REW_nobiomass

option threads=36;
option optcr= 0.05;
option MIP=gurobi;
option rmip=gurobi;
option limrow=100000;
option reslim=172800;
*option MIP=convert;
*option output=%outputname%;


sets


g                'Regions'
l                'Transportation modes'
p                'Prodution tehnologies'
r                'Reservoir'
s                'Storage tehnologies'
t                'years'
d                'Diameter size'
c                'clusters'
h                'hours'
*h                'hours'        /1*8760/
e                'renewables'
;


*Gurobi options
$onecho>gurobi.opt
         threads 36
         mipstart=1
$offecho
* mipstart=1


sets
sc(s)            'storage caverns'
sv (s)           'storage vessels'

alias (g,g1)
alias (h,h1)
;

set
N(g,g1)
Npipe(g,g1)
GR(g,r)
GS(g,s)
Gimp(g)
;



Parameters

beta                     'ratio of stored amount (%)'
y_c(p,t)                 'coefficient of CO2 capture for plant type p in time period t(tn CO2 / MWh H2)'
y_e(p,t)                 'coefficient of CO2 emission for plant type p and size j in time period t(tn CO2 / MWh H2)'
deltaH                   'Ratio of hydrogen regional pipeline operating costs to capital costs (%)'
deltaC_onshore           'ratio of onshore CO2 pipeline operating costs to capital costs'
deltaC_offshore          'ratio of offshore CO2 pipeline operating costs to capital costs'
diaH(d)                  'diameter of a regional hydrogen pipeline of diameter size d (m)'
diaC_onshore(d)          'diameter of an onshore CO2 pipeline of diameter size d (m)'
diaC_offshore(d)         'diameter of an offshore CO2 hydrogen pipeline of diameter size d (m)'
iota                     'maximum percentage of international hydrogen imports over the total demand (%)'
dur                      'duration of time periods (y)'
LTonshore                'useful life of onshore CO2 pipelines (y)'
LToffshore               'useful life of offshore CO2 pipelines (y)'
LTpipe                   'useful life of hydrogen pipelines (y)'
LTp(p)                   'useful life of hydrogen production plants (y)'
LTs(s)                   'useful life of hydrogen storage facilities (y)'
LTt(l)                   'useful life of hydrogen road transportation modes {Trailer,Tanker} (y)'
a                        'days in a year (days)'
AV(c,h,g,e)              'availability of renewable e in region g, cluster c and hour h (%)'
ayHR0(d,g,g1)            'initial availability of a regional hydrogen pipeline of diameter size d between regions g and g1 (0-1)'
ayC0(d,g,g1)             'initial availability of a onshore CO2 pipeline of diameter size d between regions g and g1 (0-1)'
aeC0(r)                  'initial availability of a offshore CO2 pipeline of diameter size d between collection point in regions g and reservoir r (0-1)'
BA(g,t)                  'biomass availability in region g and time period t'
br(g)                    'biomass regional availability'
cbio(t)                  'biomass cost in time period t (�/MWh)'
cccH(d)                  'capital costs of a regional hydrogen pipeline of diameter size q d (�k km-1)'
cccC_onshore(d)          'capital costs of an onshore CO2 pipeline of diameter size d (�k km-1)'
cccC_offshore(d)         'capital costs of an offshore CO2 pipeline of diameter size d (�k km-1)'
cgas(t)                  'natural gas cost in time period t (�/MWh)'
crf                      'capital recovery factor'
Cstart(p)                'cost for starting up for each technology type (�)'
Cshut(p)                 'cost for shutting down for each technology type (�)'
ct(t)                    'carbon tax in time period t (� kg-1 CO2)'
dc(t)                    'demand coefficient at time period t'
dem(g,t,c,h)             'total hydrogen demand in region g in time period t (MW)'
dfc(t)                   'discount factor for capital costs in time period t'
dfo(t)                   'discount factor for operating costs in time period t'
dw(l)                    'driver wage of road transportation mode l (� h-1)'
DistPipe(g,g1)           'delivery distance of a onshore CO2 pipeline between regions g and g1 (km)'
DistRes(g,r)             'Distance from CO2 collection point in region g to reservoir r (km)'
Dist(g,g1)           'regional delivery distance of hydrogen transportation mode l in region g (km)'
DistSt(g,s)              'distance between region g and underground storage type s'
DT(p)                    'min down time (h)'
ec(t)                    'cost of electricity back to grid (�/MWe)'
eta(p,t)                 'efficiency of WE in time period t (%)'
emtarget(t)              'emissions target in time period t (kgCO2)'
feR(l)                   'fuel economy of road transportation mode l transporting product type i within a region (km l-1)'
fp(l)                    'fuel price of road transportation mode l (� l-1)'
GasDem(c,h,g)            'hydrogen demand for each region g each cluster c and hour h (MWh)'
ge(l)                    'general expenses of road transportation mode l transporting product type i (� d-1)'
ir                       'discount rate (%)'
landAV(e,g)              'land availability of renewable e in region g (MW)'
lut(l)                   'load and unload time of road transportation mode l  (h)'
me(l)                    'maintenance expenses of road transportation mode l  (� km-1)'
nel                      'economic life cycle of capital investments (y)'
np0(p,g)                 'initial number of hydrogen produCtion plants of teChnology p and size j producing product type i in region g'
ns0(s,g)                 'initial number of hydrogen storage facilities of type s and size j storing product type i in region g'
pcap_max(p)              'maximum capacity of a hydrogen production plant of type p and size j (MW)'
pcap_min(p)              'manimum capacity of a hydrogen production plant of type p and size j (MW)'
pccost(p,t)              'capital cost of a production plant of type p (�k/kW)'
pimp                     'Price of hydrogen import (�/MWh)'
pocostF(p,t)             'operating production cost  in a production plant of type p  (�/MW/y)'
pocostV(p,t)             'operating production cost  in a production plant of type p  (�/MW)'
qHmax(d)                 'maximum flowrate in a hydrogen pipeline of diameter size d (kg H2 d-1)'
qCmax(d)                 'maximum flowrate in a CO2 pipeline of dimetere size d (kg H2 d-1)'
QImax(s)                 'maximum injection rate for each storage type s'
QRmax(s)                 'maximum retrieval rate for each storage type s'
rcap(r)                  'total capacity of reservoir r (kg CO2-eq)'
ri0                      'initial CO2 inventory in reservoir r (kg CO2)'
RD(p)                    'Comit Ramp down '
rccost(e,t)              'renewable e cacital cost in time period t(�/MW)'
rocost(e,t)              'renewable e operating cost in time period t(�/MW)'
RU(p)                    'Comit Ramp up '
scap_max(s)              'maximum capacity of a storage facility of type s (MWh H2)'
scap_min(s)              'minimum capacity of a storage facility of type s (MWh H2)'
sccost(s)                'capital cost of a storage facility of type s (�/MW)'
socostF(s)               'Fixed operating storage cost  in a production plant of type p (�/MW/y)'
socostV(s)               'Variable operating storage cost  in a production plant of type p (£/kWh stored)'
spR(l)                   'regional average speed of road transportation mode l transporting product type i within a region(km h-1)'
st0(s,g)                 'storage at time 0'
tcap(l)                  'capacity of road transportation mode l transporting product type i (MWh H2 unit-1)'
tmc(l)                   'capital cost of establishing a road transportation unit of transportation mode l (� unit-1)'
tmaR(l)                  'regional availability of road transportation mode l (h d-1)'
PCap(p)                  'unit capacity for production type p (MW)'
SCap(s)                  'unit capacity for storage type s (MW)'
uInit(p,g,t)             'Initial operating units type p in region g at time period t'
UT(p)                    'min up time (h)'
Vbio_max(t)              'maximum biomass consumption in year t'
WF(c)                    'weight of clusters (d)'
demH(g,t,c,h)            'demand in hour intervals (MW)'
AVH(c,h,g,e)             'availability of renewable e in region g, cluster c and hour h (%)'
;

*=============READING DATA FROM EXCEL
*$set Case 'C:\Users\Margarita\Documents\SHIPmod\newhydro'
* Read data from Excel and store in GDX file.
$call gdxxrw newhydro_clusters.xlsx skipempty=0 trace=3 index=Index!A1:F200
$gdxin newhydro_clusters.gdx
$load p,h,g,l,r,c,s,sc,sv,t,d,e,N,GR,GS,Gimp
$load beta ,y_c, y_e,dc,deltaH ,deltaC_onshore,deltaC_offshore ,diaH,diaC_onshore
$load diaC_offshore, dur , ec,emtarget, DT,LTonshore ,LToffshore ,LTpipe ,LTp, LTs, LTt
$load a,AV,ayHR0, ayC0,aeC0,br,cbio,cccH,cccC_onshore,cccC_offshore,cgas,ct,eta
$load feR,fp,GasDem,ge,landAV,lut,me,nel,np0,ns0,pcap_max,dw,DistPipe,DistRes,Dist,DistSt,PCap,pimp
$load pcap_min,pccost,pocostF,pocostV,qHmax,qCmax,QImax,QRmax,rcap,rccost,rocost,ri0,RD,RU,scap_max,SCap,scap_min,sccost
$load socostF,socostV,spR,tcap,tmc,tmaR,UT,uInit,iota,Vbio_max
$load WF
$load Cstart,Cshut
;

$gdxin



******************************************************************************************************
******************************************************************************************************
*$ONTEXT
sets
TT(t)   /3*6/
CC(c)   /1*6/
HH(h)   /1*24/;

scalar
y1       /3/
theta    /1/;
*$OFFTEXT
******************************************************************************************************
******************************************************************************************************
*HA1 SETS AND SCALARS
$ONTEXT
sets
TT(t)   /3*6/
CC(c)   /1*6/
HH(h)   /1*6/;

scalar
y1       /3/
theta    /4/;
$OFFTEXT
******************************************************************************************************
******************************************************************************************************

*electricity cost
ec(t)=3*ec(t);


ir=0.06;

*regional biomass availability
parameter bp;
bp=0.5;
display bp;
BA(g,t)=bp*br(g)*Vbio_max(t)*1000000;


*discount factor of capital costs
dfc(t)=round(1/(1+ir)**(dur*ord(t)-dur),2);
*discount factor of operating costs
dfo(t)=round(1/(1+ir)**(dur*ord(t)-5) + 1/(1+ir)**(5*ord(t)-4) + 1/(1+ir)**(5*ord(t)-3)+
1/(1+ir)**(5*ord(t)-2) + 1/(1+ir)**(5*ord(t)-1),2);
*capital recovery factor
nel=30;
crf=round((ir*(1+ir)**nel)/((1+ir)**nel-1),2);


*SET Npipe DEFINITION
Npipe(g,g1)$(DistPipe(g,g1) gt 0)=yes;
display N,Npipe;

*Set GS
st0(s,g)=0;
GS(g,sv)=yes;

*demand
dem(g,t,c,h)=round(dc(t)*GasDem(c,h,g),2);



*******************HA1*************************************************************
$ONTEXT
demH(g,t,c,h1)$HH(h1)= round(sum(h$((ord(h)<=theta*ord(h1)) and (ord(h)>theta*(ord(h1)-1))),dem(g,t,c,h)/theta),2);
AVH(c,h1,g,e)$HH(h1)=  round(sum(h$((ord(h)<=theta*ord(h1)) and (ord(h)>theta*(ord(h1)-1))),AV(c,h,g,e)/theta),2);
display demH,AVH;
$OFFTEXT


display dem;
display g,l,p,r,s,c,h,t,d,e,N,GR,GS,Gimp,beta ,y_c,y_e,Cstart, Cshut, deltaH ,deltaC_onshore,deltaC_offshore ,diaH,diaC_onshore,
diaC_offshore, dc,dur,LTonshore ,LToffshore ,LTpipe ,LTp, LTs, LTt,
ayHR0, AV,ayC0,aeC0,BA,br,cccH,cccC_onshore,cccC_offshore,cbio,cgas,crf,ct,dem,dfc,
dfo,dw,DistPipe,DistRes,Dist,DistSt,ec,emtarget,eta,feR,fp,GasDem,ge,ir,landAV,lut,me,nel,np0,ns0,pcap_max,
pcap_min,pccost,pocostF,pocostV,pimp,qHmax,qCmax,rcap,rccost,rocost,ri0,scap_max,scap_min,sccost,
socostF,socostV,spR,tcap,tmc,tmaR,pcap,scap,Vbio_max,WF;


integer Variables
InvP(p,g,t)      'Investment of new plants of type p producing in region g in time period t'
InvS(s,g,t)      'Investment of new storage facilities of type in region g in time period t'
*ITU(l,g,g1,t)    'Number of new transportation units of type l for regional transportation byroad in region g to region g acquired in time period t'
;


positive Variables
NP(p,g,t)        'Number of plants of type j and size p in region g in time period t'
NS(s,g,t)        'Number of storage facilities of type s and size p in region g in time period t'
*NTU(l,g,g1,t)    'Number of transportation units of type l for regional transportation by road in region g in time period t'
;

positive variables
ITU(l,g,g1,t)    'Number of new transportation units of type l for regional transportation byroad in region g to region g acquired in time period t'
NTU(l,g,g1,t)    'Number of transportation units of type l for regional transportation by road in region g in time period t'
;

positive Variables
AY(d,g,g1,t)     'availability of hydrogen pipelines of diameter size d for regional distribution in region g in time period t'
AYon(d,g,g1,t)   'availability of onshore CO2 pipelines of diameter size d for local distribution in region g in time period t'
AYoff(d,g,r,t)   'availability of offshore CO2 pipelines of diameter size d for local distribution in region g in time period t'
AYst(d,g,sc,t)    'availability of hydrogen pipelines of diameter size d for distribution in region g in time period t'
;

integer Variables
Yh(d,g,g1,t)     'establishment of hydrogen pipelines of diameter size d for regioanal distribution in region g in time period t'
Yon(d,g,g1,t)    'establishment of onshore CO2pipelines of diameter size d in region g in time period t'
Yoff(d,g,r,t)    'establishment of offshore CO2 pipelines of diameter size d in region g in time period t'
Yst(d,g,sc,t)     'establishment of hydrogen pipelines of diameter size d in region g in time period d to storage type s'
;

positive Variables
BC               'Biomass Cost (�)'
CL(g,t,c,h)      'curtailment in region g, time period t, cluster c, hour h (MW)'
InvR(e,g,t)      'invested capacity of renewable e in region g and time period t (MW)'
FCR              'Fuel cost for regional transport (�)'
GC               'Gas Cost (�)'
GCR              'General Cost for regional transport (�)'
IIC              'International import cost (�)'
IMP(g,t,c,h)     'flow rate of international import in region g in time period t (MW)'
LCR              'Labour cost for regional transport (�)'
MCR              'Maintenance cost for regional transport (�)'
NR(e,g,t)        'capacity of renewable e in region g and time period t (MW)'
Pr(p,g,t,c,h)    'Production rate of product type i produced by a plant of type j and size p in region g in time period t (MW)'
Pre(e,g,t,c,h)   'electricity production from renewable e in region g, time period t, cluster c, hour h (MW)'
PipeCC           'Pipeline capital cost (�k)'
PipeOC           'Pipeline operating cost (�k)'
PCC              'Production capital cost (�k)'
POC              'Production operating cost (�)'
Q(l,g,g1,t,c,h)  'regional flowrate of H2 via transportation mode l in region g in time period t (MWh)'
Qi(g,s,t,c,h)    'flowrate of H2 via pipeline from region g to storage type s in time period t(MWh)'
Qr(s,g,t,c,h)    'flowrate of H2 via pipeline from region g to storage type s in time period t(MWh)'
Qon(g,g1,t,c,h)  'regional flowrate of CO2 via onshore pipelines between regions g and g?? in time period t (kg CO2/d)'
Qoff(g,r,t,c,h)  'flowrate of CO2 via oofshore pipelines from a collection point in region g to a reservoir r in time period t (kg CO2/d)'
RCC              'Road transportation capital cost (�)'
Rdown(p,g,t,c,h) 'upward reserve contribution (MWh)'
ReC              'renewables capital cost (�)'
RI(r,t)          'Inventory of CO2 in reservoir r in time period t (kg CO2-eq)'
ROC              'Road transportation cost (�)'
Rup(p,g,t,c,h)   'downward reserve contribution (MWh)'
SCC              'Storage capical cost (�)'
SOC              'Storage operating cost (�)'
St(s,g,t,c,h)    'Average inventory of product type i stored in a storage facility of type s and size p in region g in time period t (kW)'
TCC              'Transportation capital cost (�)'
TOC              'Transportation operating cost (�)'
Vbio(t)          'Biomass consumption in time period t, cluster c and hour h (kg)'
Vgas(t)          'Gas consumpiton in timer period t, cluster c and hour h ('
slak1(g,t,c,h)
slak2(g,t,c,h)
emslak(t)
;

variables
TC               'Total cost (�)'
em(t)            'total emissions (tCO2)'
CEC              'Carbon emissions cost (�)'

     ;



*storage caverns
NS.up(s,g,t)$(sc(s) and GS(g,s))=1;
*NS.up(s,g,t)$(not(GS(g,s)))=0;


*$ontext
*Built rate for production
InvP.up('SMRCCS',g,t)$TT(t)=10;
InvP.up('ATRCCS',g,t)$TT(t)=10;
InvP.up('BECCS',g,t)$TT(t)=10;
InvP.up('WE',g,t)$TT(t)=50;
InvS.up('MPSV',g,t)$TT(t)=80;
InvS.up('HPSV',g,t)$TT(t)=80;
InvR.up(e,g,t)$TT(t)=10000;
*$offtext
*CEC.lo=-8000000000;
*pcap_min(p)=0.3*pcap_min(p);
*pcap_min(p)=0*pcap_min(p);
RU(p)$(ord(p)<=3)=0.1;
RD(p)$(ord(p)<=3)=0.1;
display  pcap_min,RU,RD;

RI.up(r,t)$TT(t)=rcap(r)/1000;
ITU.up(l,g,g1,t)$TT(t)=25;
Q.up('Pipe',g,g1,t,c,h)$(TT(t) and CC(c) and HH(h) and Npipe(g,g1))=15343;


Equations
Obj,PCCeq,POCeq,SCCeq,SOCeq,TCCeq,RCCeq,PipeCCeq,TOCeq,ROCeq,FCReq,GCReq,LCReq,MCReq,PipeOCeq,
CECeq,
IICeq,ReCeq,
h2MassBalance(g,t,c,h),co2MassBalance(g,t,c,h),
PCapacity1(p,g,t,c,h),PCapacity2(p,g,t,c,h),PAvailability(p,g,t),
SInventory2(s,g,t,c,h),MaxInj(t,g,s,c,h),MaxRetr(t,g,s,c,h),SCapacityU(s,g,t,c,h),SCapacity1(s,g,t,c,h),SCapacity2(s,g,t,c,h),SAvailability(s,g,t),
SFinal(s,g,t,c),
Transportation(g,g1,t,c,h),TAvailability(g,g1,t),
H2PipeMax(g,g1,t,c,h),OnshorePipeMax(g,g1,t,c,h),OffshorePipeMax(g,r,t,c,h),PipeStAvailability(d,g,s,t),
H2PAvailability(d,g,g1,t),OnPAvailability(d,g,g1,t),OffPAvailability(d,g,r,t),
H2Pipe(g,g1,t),OnPipe(g,g1,t),OffPipe(g,r,t),StPipe(s,g,t),
ResInventory(r,t),
ImpLimit(t,c,h),
RampUp(p,g,c,h,t),RampDown(p,g,c,h,t),
ElecProd(g,t,c,h),RenewAv(e,g,t,c,h),RenewCap(e,g,t),
CurtPercentage(c,h),
EmTargeteq(t),Emissions(t),
GasCons(t),
BioCons(t),
GasCost,BioCost,
LandAvailability(e,g,t)
BiomassAvailability(g,t)
;

*=================================OBJECTIVE FUNCTION=======================================================*
*total objective function
Obj..
TC  =e= 1000*PCC + SCC + TCC + POC + SOC + TOC + CEC + IIC + ReC + GC + BC ;
*TC  =e= 1000*PCC + SCC + TCC + POC + SOC + TOC +       IIC + ReC + GC + BC ;
   ;


*facilities capital cost
PCCeq..
PCC =e= sum((p,g,t)$TT(t),dfc(t)*pccost(p,t)*PCap(p)*InvP(p,g,t));

SCCeq..
SCC =e= sum((s,g,t)$(TT(t) and GS(g,s)),dfc(t)*sccost(s)*SCap(s)*InvS(s,g,t));

*facilities operation cost
POCeq..
POC  =e= sum((p,g,t)$TT(t),dfo(t)*(pocostF(p,t)*PCap(p)*NP(p,g,t)+ sum((c,h)$(CC(c) and HH(h)),WF(c)*(pocostV(p,t))*theta*Pr(p,g,t,c,h))));

SOCeq..
SOC  =e= sum((s,g,t)$(TT(t) and GS(g,s)),dfo(t)*(socostF(s)*SCap(s)*NS(s,g,t)+ sum((c,h)$(CC(c) and HH(h)),WF(c)*socostV(s)*theta*Qi(g,s,t,c,h))));

*sum of transportation cost
TCCeq..
TCC =e= RCC + 1000*PipeCC;

*road transportation capital cost
RCCeq..
RCC  =e=  sum((t,l,g,g1)$(TT(t) and N(g,g1) and ord(l)=1),dfc(t)*tmc(l)*ITU(l,g,g1,t));

*pipeline transportation capital cost
PipeCCeq..
PipeCC=e= sum((t,g,g1,d)$(Npipe(g,g1) and TT(t) and (ord(g)<ord(g1))),                 dfc(t)*cccH(d)*         DistPipe(g,g1)*Yh(d,g,g1,t))+
          sum((t,g,g1,d)$(N(g,g1)     and TT(t) and (ord(g)<ord(g1)) and (ord(d)<=2)), dfc(t)*cccC_onshore(d)* Dist(g,g1)*    Yon(d,g,g1,t))+
          sum((t,g,r,d)$(GR(g,r)      and TT(t) and (ord(d)<=2)),                      dfc(t)*cccC_offshore(d)*DistRes(g,r)*  Yoff(d,g,r,t))+
          sum((t,g,sc,d)$(GS(g,sc)      and TT(t)),                                    dfc(t)*cccH(d)*         DistSt(g,sc)*  Yst(d,g,sc,t));


*transportation operating cost
TOCeq..
TOC=e=ROC+ 1000*PipeOC;

*road transportation operating cost
ROCeq..
ROC =e= FCR + GCR + LCR + MCR;

*regional fuel cost
FCReq..
FCR =e= sum((t,l,g,g1,c,h)$(N(g,g1) and ord(l)=1 and TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*fp(l)*(2*Dist(g,g1)*theta*Q(l,g,g1,t,c,h)/(feR(l)*tcap(l))));

*regional general cost
GCReq..
GCR =e= sum((l,t,g,g1,c,h)$(N(g,g1) and ord(l)=1 and TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*ge(l)*theta*Q(l,g,g1,t,c,h)/tcap(l)*(2*Dist(g,g1)/spR(l)+lut(l))) ;

*regional labour cost
LCReq..
LCR =e= sum((l,t,g,g1,c,h)$(N(g,g1) and ord(l)=1 and TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*dw(l)*theta*Q(l,g,g1,t,c,h)/tcap(l)*(2*Dist(g,g1)/spR(l)+lut(l))) ;

*regional maintenance cost
MCReq..
MCR =e= sum((l,t,g,g1,c,h)$(N(g,g1) and ord(l)=1 and TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*me(l)*2*Dist(g,g1)*theta*Q(l,g,g1,t,c,h)/tcap(l)) ;

*pipline operating cost
PipeOCeq..
PipeOC=e=sum((t,l,g,g1,d)$(Npipe(g,g1) and TT(t) and (ord(g)<ord(g1))),                    dfo(t)*deltaH*         crf*cccH(d)*         DistPipe(g,g1)*AY(d,g,g1,t))
      +  sum((t,l,g,g1,d)$(N(g,g1)     and TT(t) and (ord(g)<ord(g1)) and (ord(d)<=2)),    dfo(t)*deltaC_onshore* crf*cccC_onshore(d)* Dist(g,g1)*    AYon(d,g,g1,t))
      +  sum((t,l,g,d,r)$(GR(g,r)      and TT(t) and (ord(d)<=2)) ,                        dfo(t)*deltaC_offshore*crf*cccC_offshore(d)*DistRes(g,r)*  AYoff(d,g,r,t))
      +  sum((t,g,sc,d)$(GS(g,sc) and TT(t) )  ,                       dfo(t)*deltaH*         crf*cccH(d)*         DistSt(g,sc)*  AYst(d,g,sc,t));

*carbon emissions cost
CECeq..
CEC =e= sum((p,g,t,c,h)$(TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*ct(t)*y_e(p,t)*theta*Pr(p,g,t,c,h));

*international import cost
IICeq..
IIC =e= sum((t,g,c,h)$(gimp(g)and TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*pimp*theta*IMP(g,t,c,h));

*renewables total cost
ReCeq..
ReC =E= sum((e,g,t)$TT(t), (dfc(t)*rccost(e,t)*InvR(e,g,t)+ dfo(t)*rocost(e,t)*NR(e,g,t)));

*Fuels cost
GasCost..
GC =E= sum((t)$(TT(t)), dfo(t)*cgas(t)*Vgas(t));

BioCost..
BC =E= sum((t)$(TT(t)), dfo(t)*cbio(t)*Vbio(t));

*======================================MASS BALANCES======================================================*
*hygrogen
h2MassBalance(g,t,c,h)$(TT(t) and CC(c) and HH(h))..

$ONTEXT
sum((p),Pr(p,g,t,c,h)) +
sum(g1$Npipe(g1,g), Q('Pipe',g1,g,t,c,h))+
IMP(g,t,c,h)$(gimp(g))+
sum(s$(GS(g,s)),Qr(s,g,t,c,h)) + slak2(g,t,c,h)
=E=
sum(g1$Npipe(g,g1), Q('Pipe',g,g1,t,c,h))+
sum(s$(GS(g,s)),Qi(g,s,t,c,h))+
dem(g,t,c,h) + slak1(g,t,c,h);
$OFFTEXT
******************HA2 Monolithic*************************************************
*$ONTEXT
sum((p),Pr(p,g,t,c,h)) +
sum(g1$Npipe(g1,g), Q('Pipe',g1,g,t,c,h))+
IMP(g,t,c,h)$(gimp(g))+
sum(s$(GS(g,s)),Qr(s,g,t,c,h))
=E=
sum(g1$Npipe(g,g1), Q('Pipe',g,g1,t,c,h))+
sum(s$(GS(g,s)),Qi(g,s,t,c,h))+
dem(g,t,c,h) ;
*$OFFTEXT
******************HA1*************************************************
$ONTEXT
sum((p),Pr(p,g,t,c,h)) +
sum(g1$Npipe(g1,g), Q('Pipe',g1,g,t,c,h))+
IMP(g,t,c,h)$(gimp(g))+
sum(s$(GS(g,s)),Qr(s,g,t,c,h))
=E=
sum(g1$Npipe(g,g1), Q('Pipe',g,g1,t,c,h))+
sum(s$(GS(g,s)),Qi(g,s,t,c,h))+
demH(g,t,c,h) ;
$OFFTEXT


*sum((p),Pr(p,g,t,c,h)) +
*sum((g1)$(N(g1,g)),Q('Trailer',g1,g,t,c,h)) +
*sum(g1$Npipe(g1,g), Q('Pipe',g1,g,t,c,h))+
*IMP(g,t,c,h)$(gimp(g))+
*sum(s$(GS(g,s)),Qr(s,g,t,c,h))
*=e=
*sum((g1)$(N(g,g1)),Q('Trailer',g,g1,t,c,h))  +
*sum(g1$Npipe(g,g1), Q('Pipe',g,g1,t,c,h))+
*sum(s$(GS(g,s)),Qi(g,s,t,c,h))+
*dem(g,t,c,h);

*co2
co2MassBalance(g,t,c,h)$(TT(t) and CC(c) and HH(h))..

sum(g1$N(g1,g),Qon(g1,g,t,c,h))+sum(p,y_c(p,t)*Pr(p,g,t,c,h)) =e=
sum(g1$N(g,g1),Qon(g,g1,t,c,h))+sum(r$GR(g,r),Qoff(g,r,t,c,h));

*======================================FUELS CONSUMPTION========================================================*
*Gas consumption
GasCons(t)$(TT(t))..
Vgas(t) =e= sum((g,p,c,h)$(ord(p)<=2 and CC(c) and HH(h)),WF(c)*theta*Pr(p,g,t,c,h)/eta(p,t)) ;

*Biomass consumption
BioCons(t)$(TT(t))..
Vbio(t) =e= sum((g,c,h)$(CC(c) and HH(h)),WF(c)*theta*Pr('BECCS',g,t,c,h)/eta('BECCS',t)) ;

*Biomass availability
BiomassAvailability(g,t)$(TT(t))..
sum((c,h)$(CC(c) and HH(h)),WF(c)*theta*Pr('BECCS',g,t,c,h)/eta('BECCS',t)) =l= BA(g,t);
*======================================RAMP UP/DOWN========================================================*

*Ramp Up
RampUp(p,g,c,h,t)$(TT(t) and CC(c) and HH(h) and ord(h)>1)..
Pr(p,g,t,c,h) - Pr(p,g,t,c,h-1)  =l= theta*RU(p)*PCap(p)*NP(p,g,t);


*Ramp Down
RampDown(p,g,c,h,t)$(TT(t) and CC(c) and HH(h) and ord(h)>1)..
Pr(p,g,t,c,h-1)-Pr(p,g,t,c,h)  =l= theta*RD(p)*PCap(p)*NP(p,g,t);

*======================================PRODUCTION CONSTRAINTS==============================================*
PCapacity1(p,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr(p,g,t,c,h)=g=PCap(p)*pcap_min(p)*NP(p,g,t)  ;


PCapacity2(p,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr(p,g,t,c,h)=l=PCap(p)*pcap_max(p)*NP(p,g,t)  ;


*availability in production plants
PAvailability(p,g,t)$(TT(t))..
*NP(p,g,t) =e= NP(p,g,t-1)$(ord(t)>y1) + np0(p,g)$(ord(t)=y1) + InvP(p,g,t)  - InvP(p,g,t-(LTp(p)/dur))$((ord(t)-(LTp(p)/dur)) gt y1);
NP(p,g,t) =e= NP(p,g,t-1)$(ord(t)>y1) + np0(p,g)$(ord(t)=y1) + InvP(p,g,t) ;

*======================================STORAGE CONSTRAINTS=================================================*
*total average inventory for underground storage
SInventory2(s,g,t,c,h)$(GS(g,s)and TT(t) and CC(c) and HH(h))..
St(s,g,t,c,h) =e= St(s,g,t,c,h-1)$(ord(h)>1) + st0(s,g)$(ord(h)=1) +theta*(Qi(g,s,t,c,h) - Qr(s,g,t,c,h)) ;

*maximum injection and retrieval rate
MaxInj(t,g,s,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
Qi(g,s,t,c,h)=l=QImax(s)*NS(s,g,t);

MaxRetr(t,g,s,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
Qr(s,g,t,c,h)=l=QRmax(s)*NS(s,g,t);

*storage capacity constraints for underground storage
SCapacityU(sc,g,t,c,h)$(GS(g,sc) and TT(t) and CC(c) and HH(h))..
InvS(sc,g,t)=l=sum(d,Yst(d,g,sc,t));

*storage capacity constraints
SCapacity1(s,g,t,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
St(s,g,t,c,h)=g=SCap(s)*scap_min(s)*NS(s,g,t)  ;

SCapacity2(s,g,t,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
St(s,g,t,c,h)=l=SCap(s)*scap_max(s)*NS(s,g,t) ;

*availability of storage facility
SAvailability(s,g,t)$(GS(g,s) and TT(t))..
NS(s,g,t) =e=  NS(s,g,t-1)$(ord(t)>y1) +  ns0(s,g)$(ord(t)=y1) +InvS(s,g,t);

SFinal(s,g,t,c)$(GS(g,s) and TT(t) and CC(c))..
St(s,g,t,c,'24') =e= 0;
*======================================RENEWABLES CONSTRAINTS==============================================*
*Electricity production for electrolysis
ElecProd(g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr('WE',g,t,c,h) =e= eta('WE',t)*(sum(e,Pre(e,g,t,c,h)) - CL(g,t,c,h));

*Renewables Availability
RenewAv(e,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pre(e,g,t,c,h) =e= 0.7*AV(c,h,g,e)*NR(e,g,t);
******************HA1*************************************************
*Pre(e,g,t,c,h) =e= 0.7*AVH(c,h,g,e)*NR(e,g,t);

*Renewables capacity expansion
RenewCap(e,g,t)$(TT(t))..
NR(e,g,t) =e= NR(e,g,t-1)$(ord(t)>y1) + InvR(e,g,t);


*Land availability
LandAvailability(e,g,t)$(TT(t))..
NR(e,g,t) =l= landAV(e,g);

*Curtailment limit
CurtPercentage(c,h)$(CC(c) and HH(h))..
sum((g,t),CL(g,t,c,h)) =l= 0.1*sum((e,g,t), Pre(e,g,t,c,h)) ;

*=====================================ROAD TRANSPORTATION CONSTRAINTS=========================================*
*regional transportation
Transportation(g,g1,t,c,h)$(N(g,g1) and TT(t) and CC(c) and HH(h))..
NTU('Trailer',g,g1,t) =g= Q('Trailer',g,g1,t,c,h)/(tcap('Trailer'));

*availability of transportation
TAvailability(g,g1,t)$(N(g,g1) and 1 and TT(t))..
NTU('Trailer',g,g1,t) =e=  NTU('Trailer',g,g1,t-1)$(ord(t)>y1)+ ITU('Trailer',g,g1,t);

*=====================================PIPELINE CONSTRAINTS=====================================================*
*maximum flowrate for pipelines
*H2PipeMax(g,g1,t,c,h)$(Npipe(g,g1) and ord(t)<=tm and ord(c)<=cm and ord(h)<=hm)..      Q('Pipe',g,g1,t,c,h)=l= sum(d,qHmax(d)*AY(d,g,g1,t));
*OnshorePipeMax(g,g1,t,c,h)$(N(g,g1) and ord(t)<=tm and ord(c)<=cm and ord(h)<=hm)..     Qon(g,g1,t,c,h)=l= sum(d,qCmax(d)*AYon(d,g,g1,t));
H2PipeMax(g,g1,t,c,h)$(Npipe(g,g1) and TT(t) and CC(c) and HH(h))..
Q('Pipe',g,g1,t,c,h)=l= sum(d,qHmax(d)*(AY(d,g,g1,t)$(ord(g)<ord(g1))+AY(d,g1,g,t)$(ord(g1)<ord(g))));

OnshorePipeMax(g,g1,t,c,h)$(N(g,g1) and TT(t) and CC(c) and HH(h))..
*Qon(g,g1,t,c,h)=l= sum(d,10*qCmax(d)*(AYon(d,g,g1,t)$(ord(g)<ord(g1))+AYon(d,g1,g,t)$(ord(g1)<ord(g))));
Qon(g,g1,t,c,h)=l= sum(d$(ord(d)<=2),qCmax(d)*(AYon(d,g,g1,t)$(ord(g)<ord(g1))+AYon(d,g1,g,t)$(ord(g1)<ord(g))));

OffshorePipeMax(g,r,t,c,h)$(GR(g,r) and TT(t) and CC(c) and HH(h))..
*Qoff(g,r,t,c,h)=l= sum(d,10*qCmax(d)*AYoff(d,g,r,t));
Qoff(g,r,t,c,h)=l= sum(d$(ord(d)<=2),qCmax(d)*AYoff(d,g,r,t));

*availability of pipelines
H2PAvailability(d,g,g1,t)$(Npipe(g,g1) and (ord(g)<ord(g1)) and TT(t))..
*AY(d,g,g1,t)   =e= AY(d,g,g1,t-1)$(ord(t)>y1)  +  ayHR0(d,g,g1)$(ord(t)=y1) + Yh(d,g,g1,t) - Yh(d,g,g1,t-(LTpipe/dur))$(ord(t)-(LTpipe/dur) gt y1)  ;
AY(d,g,g1,t)   =e= AY(d,g,g1,t-1)$(ord(t)>y1)  +  ayHR0(d,g,g1)$(ord(t)=y1) + Yh(d,g,g1,t) ;

OnPAvailability(d,g,g1,t)$(N(g,g1) and (ord(g)<ord(g1)) and TT(t) and ord(d)<=2)..
*AYon(d,g,g1,t) =e= AYon(d,g,g1,t-1)$(ord(t)>y1) +  ayC0(d,g,g1)$(ord(t)=y1) + Yon(d,g,g1,t) - Yon(d,g,g1,t-(LTonshore/dur))$((ord(t)-(LTonshore/dur)) gt y1)  ;
AYon(d,g,g1,t) =e= AYon(d,g,g1,t-1)$(ord(t)>y1) +  ayC0(d,g,g1)$(ord(t)=y1) + Yon(d,g,g1,t) ;

OffPAvailability(d,g,r,t)$(GR(g,r) and TT(t) and ord(d)<=2)..
*AYoff(d,g,r,t) =e= AYoff(d,g,r,t-1)$(ord(t)>y1) +  aeC0(r)$(ord(t)=1) + Yoff(d,g,r,t) - Yoff(d,g,r,t-(LToffshore/dur))$((ord(t)-(LToffshore/dur)) gt y1)   ;
AYoff(d,g,r,t) =e= AYoff(d,g,r,t-1)$(ord(t)>y1) +  aeC0(r)$(ord(t)=y1) + Yoff(d,g,r,t) ;

PipeStAvailability(d,g,sc,t)$(GS(g,sc) and TT(t))..
*AYst(d,g,sc,t) =e= AYst(d,g,sc,t-1)$(ord(t)>y1) +  0$(ord(t)=1) + Yst(d,g,sc,t) - Yst(d,g,sc,t-(LTpipe/dur))$((ord(t)-(LTpipe/dur)) gt y1)   ;
AYst(d,g,sc,t) =e= AYst(d,g,sc,t-1)$(ord(t)>y1) + Yst(d,g,sc,t) ;

*one diameter size
H2Pipe(g,g1,t)$(Npipe(g,g1) and (ord(g)<ord(g1)) and TT(t))..
sum(d,AY(d,g,g1,t))=l=1;

OnPipe(g,g1,t)$(N(g,g1) and (ord(g)<ord(g1)) and TT(t))..
sum(d$(ord(d)<=2),AYon(d,g,g1,t))=l=1;

OffPipe(g,r,t)$(GR(g,r)and TT(t))..
sum(d$(ord(d)<=2),AYoff(d,g,r,t))=l=1;

StPipe(sc,g,t)$(GS(g,sc)and TT(t))..
sum(d,AYst(d,g,sc,t))=l=1;


*======================================RESERVOIRS==============================================================*
*inventory
ResInventory(r,t)$TT(t)..
RI(r,t)=e=RI(r,t-1)$(ord(t)>=y1) +ri0(r)$(ord(t)=y1)/1000 + dur*sum((g,c,h)$GR(g,r),WF(c)*theta*Qoff(g,r,t,c,h))$(ord(t)>=y1)/1000;

*inventory level maximum
*ResCapacity(r,t)$TT(t)..                      RI(r,t)=l=sum((d,g)$GR(g,r),AYoff(d,g,r,t)*rcap(r));


*======================================IMPORT=================================================================*
*Import limit
ImpLimit(t,c,h)$(TT(t) and CC(c) and HH(h))..
sum(g$Gimp(g),Imp(g,t,c,h))=l=iota*sum(g,dem(g,t,c,h));


*======================================EMISSIONS================================================================*
*Emissions target
Emissions(t)$TT(t)..
em(t) =e= sum((p,g,c,h)$(CC(c) and HH(h)),WF(c)*y_e(p,t)*theta*Pr(p,g,t,c,h));

EmTargeteq(t)$TT(t)..
em(t) =l= emtarget(t) ;


*======================================TIGHTENING===============================================================*
Equations tightGas,tightBio, tightinvWE, tightInvStorage, tightITU;

tightGas(t)$( TT(t))..
sum((g,p)$(ord(p) <= 2), INVP(p,g,t)) =l= 40;

tightBio(p,t)$(ord(p) = 3 and TT(t))..
sum(g, INVP(p,g,t)) =l= 30;

tightInvWE(p,t)$(ord(p) = card(p) and TT(t))..
sum(g, INVP(p,g,t)) =l= 200;

tightInvStorage(s,t)$(ord(s)>4 and TT(t))..
sum(g, INVS(s,g,t)) =l= 5*50;



*tightITU(t,g1)..
*sum(g, ITU('Trailer',g,g1,t)) =l= 250;

******************************BOUNDS********************************************
*Vgas.up(t)$TT(t)=50000000;
*https://www.tolvik.com/wp-content/uploads/2020/04/Tolvik-UK-Biomass-Statistics-2019-FINAL.pdf
*Vbio.up(t)$TT(t)=50000000;

$ontext
Vbio.up('1')=21210000;
Vbio.up('1')=25370000;
Vbio.up('1')=29530000;
Vbio.up('1')=33690000;
Vbio.up('1')=37840000;
Vbio.up('1')=42000000;
$offtext


Qon.up(g,g1,t,c,h)$(TT(t) and CC(c) and HH(h))=1.17E+04 ;
Qoff.up(g,r,t,c,h)$(TT(t) and CC(c) and HH(h))=1.17E+04 ;
*Pr.up(p,g,t,c,h)$(TT(t) and CC(c) and HH(h))=20000;
*St.up(s,g,t,c,h)$(TT(t) and CC(c) and HH(h))=80000;

model hydrogen /Obj,PCCeq,POCeq,SCCeq,SOCeq,TCCeq,PipeCCeq,TOCeq,PipeOCeq,
CECeq,
IICeq,ReCeq,
*RCCeq,ROCeq,FCReq,GCReq,LCReq,MCReq,
h2MassBalance,co2MassBalance,
*PCapacity1,
PCapacity2,PAvailability,
SInventory2,MaxInj,MaxRetr,SCapacityU,SCapacity1,SCapacity2,SAvailability,
SFinal,
*Transportation,TAvailability,
H2PipeMax,OnshorePipeMax,OffshorePipeMax,PipeStAvailability,
H2PAvailability,OnPAvailability,OffPAvailability,
H2Pipe,OnPipe,OffPipe,StPipe,
ResInventory,
ImpLimit,
RampUp,RampDown,
ElecProd,RenewAv,RenewCap,
*CurtPercentage,
*tightGas, tightBio,tightinvWE,
*tightInvStorage
EmTargeteq,
Emissions,
GasCons,
BioCons,
GasCost,BioCost,
LandAvailability,
BiomassAvailability
/
;

model hydrogenHA2 /Obj,PCCeq,POCeq,SCCeq,SOCeq,TCCeq,PipeCCeq,TOCeq,PipeOCeq,
CECeq,
IICeq,ReCeq,
*RCCeq,ROCeq,FCReq,GCReq,LCReq,MCReq,
h2MassBalance,co2MassBalance,
*PCapacity1,
PCapacity2,PAvailability,
SInventory2,MaxInj,MaxRetr,SCapacityU,SCapacity1,SCapacity2,SAvailability,
SFinal,
*Transportation,TAvailability,
*H2PipeMax,OnshorePipeMax,OffshorePipeMax,
H2PAvailability,OnPAvailability,OffPAvailability,PipeStAvailability,
H2Pipe,OnPipe,OffPipe,StPipe,
ResInventory,
ImpLimit,
RampUp,RampDown,
ElecProd,RenewAv,RenewCap,
*CurtPercentage,
*tightGas, tightBio, tightinvWE,
*tightInvStorage
EmTargeteq,
Emissions,
GasCons,
BioCons,
GasCost,
BioCost,
LandAvailability,
BiomassAvailability
/
;

hydrogenHA2.optfile=1;
hydrogen.optfile=1;



*scaling
PCCeq.scale=1000;
POCeq.scale=1000;
SCCeq.scale=1000;
SOCeq.scale=1000;
ReCeq.scale=1000;
PipeCCeq.scale=1000;
PipeOCeq.scale=1000;
OnshorePipeMax.scale(g,g1,t,c,h)=1000;
OffshorePipeMax.scale(g,r,t,c,h)=1000;
EmTargeteq.scale(t)=10000;
BiomassAvailability.scale(g,t)=1000;
LandAvailability.scale(e,g,t)=1000;
SCapacity2.scale(s,g,t,c,h)=1000;

hydrogen.scaleopt=1;
hydrogenha2.scaleopt=1;

********************************************************************************
***********************************HA2******************************************
*$ontext
solve hydrogenHA2 using MIP minimizing TC ;


display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,RCC.l,RI.l,St.l,
*CEC.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC.l,ROC.l,SCC.l,SOC.l,TCC.l,TOC.l,
Vbio.l,Vgas.l,
BC.l,GC.l,TC.l;
*slak1.l,slak2.l,emslak.l
*ITU.l,NTU.l,RCC.l,ROC.l,FCR.l,GCR.l,LCR.l,MCR.l;


InvP.fx(p,g,t)$TT(t)=round(InvP.l(p,g,t));
InvS.fx(s,g,t)$TT(t)=round(InvS.l(s,g,t));
*ITU.fx(l,g,g1,t)$TT(t)=fixITU(l,g,g1,t);
*Yh.fx(d,g,g1,t)$TT(t)=fixYh(d,g,g1,t);
*Yon.fx(d,g,g1,t)$TT(t)=fixYon(d,g,g1,t);
*Yoff.fx(d,g,r,t)$TT(t)=fixYoff(d,g,r,t);
*Yst.fx(d,g,s,t)$TT(t)=fixYst(d,g,s,t);



*$offtext

********************************************************************************
***********************************HA1******************************************
$ontext
solve hydrogen using MIP minimizing TC ;


display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,RCC.l,RI.l,St.l,
Vbio.l,Vgas.l,
*CEC.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC.l,ROC.l,SCC.l,SOC.l,TCC.l,TOC.l,
BC.l,GC.l,TC.l;

*ITU.l,NTU.l,RCC.l,ROC.l,FCR.l,GCR.l,LCR.l,MCR.l;


InvP.fx(p,g,t)=round(InvP.l(p,g,t));
InvS.fx(s,g,t)=round(InvS.l(s,g,t));
*ITU.fx(l,g,g1,t)=round(ITU(l,g,g1,t));
*Yh.fx(d,g,g1,t)=round(Yh.l(d,g,g1,t));
*Yon.fx(d,g,g1,t)=round(Yon.l(d,g,g1,t));
*Yoff.fx(d,g,r,t)=round(Yoff.l(d,g,r,t));
*Yst.fx(d,g,sc,t)=round(Yst.l(d,g,sc,t));


HH(h)$(ord(h)<=24)=yes;
theta=1;
demH(g,t,c,h1)$HH(h1)= round(sum(h$((ord(h)<=theta*ord(h1)) and (ord(h)>theta*(ord(h1)-1))),dem(g,t,c,h)),2);
AVH(c,h1,g,e)$HH(h1)= round(sum(h$((ord(h)<=theta*ord(h1)) and (ord(h)>theta*(ord(h1)-1))),AV(c,h1,g,e)),2);

$offtext

********************************************************************************
*******************************WARM START***************************************
$ontext
solve hydrogenHA2 using MIP minimizing TC ;


display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,RCC.l,RI.l,St.l,
*CEC.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC.l,ROC.l,SCC.l,SOC.l,TCC.l,TOC.l,
Vbio.l,Vgas.l,
BC.l,GC.l,TC.l;
*slak1.l,slak2.l,emslak.l
*ITU.l,NTU.l,RCC.l,ROC.l,FCR.l,GCR.l,LCR.l,MCR.l;


InvP.fx(p,g,t)$TT(t)=round(InvP.l(p,g,t));
InvS.fx(s,g,t)$TT(t)=round(InvS.l(s,g,t));
*ITU.fx(l,g,g1,t)$TT(t)=fixITU(l,g,g1,t);
*Yh.fx(d,g,g1,t)$TT(t)=fixYh(d,g,g1,t);
*Yon.fx(d,g,g1,t)$TT(t)=fixYon(d,g,g1,t);
*Yoff.fx(d,g,r,t)$TT(t)=fixYoff(d,g,r,t);
*Yst.fx(d,g,s,t)$TT(t)=fixYst(d,g,s,t);

solve hydrogen using MIP minimizing TC ;

display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,RCC.l,RI.l,St.l,
*CEC.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC.l,ROC.l,SCC.l,SOC.l,TCC.l,TOC.l,
Vbio.l,Vgas.l,
BC.l,GC.l,TC.l;
*slak1.l,slak2.l,emslak.l
*ITU.l,NTU.l,RCC.l,ROC.l,FCR.l,GCR.l,LCR.l,MCR.l;

*parameters InvPpar(p,g,t),InvSpar(s,g,t),Yhpar(d,g,g1,t),Yonpar(d,g,g1,t),
*Yoffpar(d,g,r,t),Ystpar(d,g,sc,t);


*InvPpar(p,g,t)$TT(t)=round(InvP.l(p,g,t));
*InvSpar(s,g,t)$TT(t)=round(InvS.l(s,g,t));
*Yhpar(d,g,g1,t)$TT(t)=round(Yh.l(d,g,g1,t));
*Yonpar(d,g,g1,t)$TT(t)=round(Yon.l(d,g,g1,t));
*Yoffpar(d,g,r,t)$TT(t)=round(Yoff.l(d,g,r,t));
*Ystpar(d,g,sc,t)$TT(t)=round(Yst.l(d,g,sc,t));

*InvP.l(p,g,t)$TT(t)=InvPpar(p,g,t);
*InvS.l(s,g,t)$TT(t)=InvSpar(s,g,t);
*Yh.l(d,g,g1,t)$TT(t)=Yhpar(d,g,g1,t);
*Yon.l(d,g,g1,t)$TT(t)=Yonpar(d,g,g1,t);
*Yoff.l(d,g,r,t)$TT(t)=Yoffpar(d,g,r,t);
*Yst.l(d,g,sc,t)$TT(t)=Ystpar(d,g,sc,t);


$offtext

********************************************************************************
******************************WARM START II*************************************
$ontext
CC(c)=no;
CC('1')=yes;
CC('2')=yes;
display CC;

option optcr= 0.1;

solve hydrogen using MIP minimizing TC ;

display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,RCC.l,RI.l,St.l,
*CEC.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC.l,ROC.l,SCC.l,SOC.l,TCC.l,TOC.l,
Vbio.l,Vgas.l,
BC.l,GC.l,TC.l;
*slak1.l,slak2.l,emslak.l
*ITU.l,NTU.l,RCC.l,ROC.l,FCR.l,GCR.l,LCR.l,MCR.l;

parameters InvPpar(p,g,t),InvSpar(s,g,t),Yhpar(d,g,g1,t),Yonpar(d,g,g1,t),
Yoffpar(d,g,r,t),Ystpar(d,g,sc,t);


InvPpar(p,g,t)$TT(t)=round(InvP.l(p,g,t));
InvSpar(s,g,t)$TT(t)=round(InvS.l(s,g,t));
Yhpar(d,g,g1,t)$TT(t)=round(Yh.l(d,g,g1,t));
Yonpar(d,g,g1,t)$TT(t)=round(Yon.l(d,g,g1,t));
Yoffpar(d,g,r,t)$TT(t)=round(Yoff.l(d,g,r,t));
Ystpar(d,g,sc,t)$TT(t)=round(Yst.l(d,g,sc,t));

InvP.l(p,g,t)$TT(t)=InvPpar(p,g,t);
InvS.l(s,g,t)$TT(t)=InvSpar(s,g,t);
*Yh.l(d,g,g1,t)$TT(t)=Yhpar(d,g,g1,t);
*Yon.l(d,g,g1,t)$TT(t)=Yonpar(d,g,g1,t);
*Yoff.l(d,g,r,t)$TT(t)=Yoffpar(d,g,r,t);
*Yst.l(d,g,sc,t)$TT(t)=Ystpar(d,g,sc,t);

CC(c)$(ord(c)<=6)=yes;
display CC;

option optcr= 0.05;

$offtext

********************************************************************************
*******************************MONOLITHIC***************************************
option optcr= 0.05;
solve hydrogen using MIP minimizing TC ;

display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,RCC.l,RI.l,St.l,
Vbio.l,Vgas.l,
CEC.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC.l,ROC.l,SCC.l,SOC.l,TCC.l,TOC.l,
BC.l,GC.l,TC.l;
*slak1.l,slak2.l,emslak.l


$ONTEXT
***************************************************************************************
*****24H MODEL RUN*******************************************************************



parameters
HourlyGas(h,g)
HourlyAV(h,g,e);

$call gdxxrw HourlyHeatDemand.xlsx skipempty=0 trace=3 index=Index!A1:F200
$gdxin HourlyHeatDemand.gdx
$load HourlyGas,HourlyAV
;
$gdxin


HH(h)$(ord(h)<=8760)=yes;
CC(c)=no;
CC(c)$(ord(c)=1)=yes;

dem(g,t,c,h)=0;
dem(g,t,c,h)$(CC(c) and HH(h)) = round(dc(t)*HourlyGas(h,g),2);


AV(c,h,g,e)=0;
AV(c,h,g,e)$(CC(c) and HH(h)) =  HourlyAV(h,g,e);

display  HourlyGas,HourlyAV,dem,AV;

InvP.fx(p,g,t)$TT(t)=round(InvP.l(p,g,t));
InvS.fx(s,g,t)$TT(t)=round(InvS.l(s,g,t));
*ITU.fx(l,g,g1,t)$TT(t)=round(ITU.l(l,g,g1,t));
Yh.fx(d,g,g1,t)$TT(t)=round(Yh.l(d,g,g1,t));
Yon.fx(d,g,g1,t)$TT(t)=round(Yon.l(d,g,g1,t));
Yoff.fx(d,g,r,t)$TT(t)=round(Yoff.l(d,g,r,t));
Yst.fx(d,g,sc,t)$TT(t)=round(Yst.l(d,g,sc,t));

solve hydrogen using MIP minimizing TC ;

display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
CEC.l,IIC.l,IMP.l,Pr.l,Pre.l,PipeCC.l,PipeOC.l,PCC.l,POC.l,Q.l,Qi.l,Qr.l,Qon.l,
Qoff.l,RCC.l,RI.l,ROC.l,SCC.l,SOC.l,St.l,TCC.l,TOC.l,TC.l;

$OFFTEXT

********************************************************************************
********************************************************************************



parameters prodcapacity(p,g,TT),prodcapacityperregion(g,TT),prodcapacitypertech(p,TT), storcapacity(s,g,TT),storcapacityperregion(g,TT);
prodcapacity(p,g,TT)=NP.l(p,g,TT)*Pcap(p) + eps$[ NOT NP.l(p,g,TT)];
prodcapacityperregion(g,TT)=sum(p,prodcapacity(p,g,TT));
prodcapacitypertech(p,TT)= sum(g,prodcapacity(p,g,TT));
storcapacity(s,g,TT)=NS.l(s,g,TT)*Scap(s) + eps$[ NOT NS.l(s,g,TT)];
storcapacityperregion(g,TT)=sum(s,storcapacity(s,g,TT));

parameters prodprofil(p,g,TT,c,h),storprofil(s,g,TT,c,h);
prodprofil(p,g,TT,CC,HH) = Pr.l(p,g,TT,CC,HH)+ eps$[ NOT Pr.l(p,g,TT,CC,HH)];
storprofil(s,g,TT,CC,HH)$GS(g,s) =St.l(s,g,TT,CC,HH)+ eps$[ NOT St.l(s,g,TT,CC,HH)];

parameters H2flow(l,g,g1,TT,c,h),Co2onflow(g,g1,TT,c,h),Co2offflow(g,r,TT,c,h);
H2flow(l,g,g1,TT,CC,HH)=  Q.l(l,g,g1,TT,CC,HH) + eps$[ NOT Q.l(l,g,g1,TT,CC,HH)];
Co2onflow(g,g1,TT,CC,HH)= Qon.l(g,g1,TT,CC,HH) + eps$[ NOT Qon.l(g,g1,TT,CC,HH)];
Co2offflow(g,r,TT,CC,HH)= Qoff.l(g,r,TT,CC,HH) + eps$[ NOT Qoff.l(g,r,TT,CC,HH)];

display prodcapacity,prodcapacityperregion,storcapacity,storcapacityperregion,prodprofil,storprofil;

parameter H2cost(TT),levelised;
H2cost(TT) = sum((g,CC,HH),WF(CC)*h2MassBalance.m(g,TT,CC,HH))/(card(g)*card(CC)*card(HH));
levelised=TC.l/(sum((g,TT,CC,HH),dfo(TT)*WF(CC)*dem(g,TT,CC,HH)));
display H2cost,levelised,h2MassBalance.m;

parameter totalimport(TT),per_import(TT);
totalimport(TT)=sum((g,CC,HH),WF(CC)*IMP.l(g,TT,CC,HH));
per_import(TT)=totalimport(TT)/sum((g,CC,HH),WF(CC)*dem(g,TT,CC,HH)) ;
display totalimport,per_import;


parameter totdem(TT,c,h),totprod(TT,c,h),totstor(TT,c,h);
totdem(TT,CC,HH)=sum((g),dem(g,TT,CC,HH))/1000;
totprod(TT,CC,HH)=sum((p,g),Pr.l(p,g,TT,CC,HH))/1000 ;
totstor(TT,CC,HH)=sum((s,g)$GS(g,s),St.l(s,g,TT,CC,HH))/1000 ;
display totdem,totprod,totstor;

parameter Vbioregional(g,TT);
Vbioregional(g,TT)= sum((c,h)$(CC(c) and HH(h)),WF(c)*theta*Pr.l('BECCS',g,TT,c,h)/eta('BECCS',TT));
display Vbioregional;

parameter techload(p,TT), techloadregion(p,g,TT),techloadcluster(p,c,TT);
techload(p,TT)= sum((g,CC,HH),WF(CC)*Pr.l(p,g,TT,CC,HH))/(PCap(p)*sum(g,NP.l(p,g,TT))*8760);
techloadregion(p,g,TT)= sum((CC,HH),WF(CC)*Pr.l(p,g,TT,CC,HH))/(PCap(p)*NP.l(p,g,TT)*8760);
techloadcluster(p,c,TT)$CC(c)= sum((g,HH),Pr.l(p,g,TT,c,HH))/(24*PCap(p)*sum(g,NP.l(p,g,TT)));
parameter techload, techloadregion,techloadcluster;

parameter
param(*)
;

param("Production Capital Cost") = 1000*PCC.l/1E9;
param("Production Operating Cost") = POC.l/1E9;
param("Storage Capital Cost") = SCC.l/1E9;
param("Storage Operating Cost") = SOC.l/1E9;
param("Pipeline Cost") = 1000*(PipeCC.l + PipeOC.l)/1E9;
*param("Road Transportation Cost") = RCC.l + ROC.l;
param("Carbon Emissions Cost") = CEC.l/1E9;
param("Imports Cost") = IIC.l/1E9;
param("Renewables Cost")=ReC.l/1E9;
param("Biomass Cost")=BC.l/1E9;
param("Gas Cost")=GC.l/1E9;
param("Total Cost")=TC.l/1E9;

*EXTRACT TO EXCEL
execute_unload  "%outputname%.gdx",
prodcapacity,
prodcapacityperregion,
prodcapacitypertech,
storcapacity,
storcapacityperregion,
prodprofil,
storprofil,
H2flow,
Co2onflow,
Co2offflow,
AY.l,
AYon.l,
AYoff.l,
Q.l,Qon.l,Qoff.l
em.l,
emtarget,
param,
H2cost,
levelised,
dem,
totalimport,
per_import
totdem,totprod,totstor,
Vbioregional,BA,
techload, techloadregion,techloadcluster;


$onecho >  exceloutput.txt

epsout=0
par=prodcapacity rng=prodcapacity!                       Rdim=2 Cdim=1
par=prodcapacityperregion rng=prodcapacityperregion!     Rdim=1 Cdim=1
par=prodcapacitypertech rng=prodcapacitypertech!         Rdim=1 Cdim=1
par=storcapacity rng=storcapacity!                       Rdim=2 Cdim=1
par=storcapacityperregion rng=storcapacityperregion!     Rdim=1 Cdim=1
par=prodprofil rng=prodprofil!                           Rdim=3 Cdim=2
par=storprofil rng=storprofil!                           Rdim=3 Cdim=2
var=AY.L rng=Pipelines!A1                                Rdim=3 Cdim=1
var=AYon.L rng=Pipelines!I1                              Rdim=3 Cdim=1
var=AYoff.L rng=Pipelines!Q1                             Rdim=3 Cdim=1
par=H2flow rng=H2flows!                                  Rdim=5 cdim=1
par=Co2onflow rng=CO2onflows!                            Rdim=4 cdim=1
par=Co2offflow rng=CO2offflows!                          Rdim=4 cdim=1
par=param rng=Costs!                                     Rdim=1
par=levelised rng=H2levelised!
par=H2cost  rng=H2levelised!A4                           Rdim=1
var=em.l    rng=Emissions!A1                             Cdim=1
par=emtarget    rng=Emissions!H1                         Cdim=1
par=dem     rng=Demand!                                  Rdim=3 cdim=1
par=totalimport rng=imports!A1                           Cdim=1
par=per_import rng=imports!A5                            Cdim=1
par=totdem       rng=totprofiles!A1                      rdim=2 cdim=1
par=totprod      rng=totprofiles!A40                     rdim=2 cdim=1
par=totstor      rng=totprofiles!a80                     rdim=2 cdim=1
par=Vbioregional rng=Vbioregional!a1                     rdim=1 cdim=1
par=BA           rng=Vbioregional!a18                    rdim=1 cdim=1
par=techload     rng=techload!A1                         rdim=1 cdim=1
par=techloadregion   rng=techload!A10                    rdim=1 cdim=2
par=techloadcluster      rng=techload!A20                rdim=1 cdim=2
$offecho

Execute 'gdxxrw %outputname%.gdx output=%outputname%.xlsx @exceloutput';

*$offtext
