#
# Plot distance difference matrix
##################################################################
##########################SETTINGS################################
moleculeA='C:/Users/brain/Downloads/Re__phyA_XFEL_update/DDM/4511_Pr_filtered.pdb'
moleculeB='C:/Users/brain/Downloads/Re__phyA_XFEL_update/DDM/4511_Pfr_filtered.pdb'
chainA='B'
chainB='B'
maxdifference=1.0

# The program will plot a distance difference matrix for Calpha atoms
# (Molecule B) - (Molecule A)
# The scale of the color map runs from -maxdifference to +maxdifference

##################################################################
##################################################################
##################################################################
##################################################################


# import modules
import string
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# set matplotlib parameter such that pdf files can be used in
# adobe illustrator etc.
matplotlib.rcParams['pdf.fonttype'] = 42
font = {'family' : 'sans-serif', 'weight' : 'normal', 'sans-serif' : "Arial", "size":10  } # 'size' : 7
matplotlib.rc('font', **font) 




################Read a pdb file###################################
def readpdb(infile,chain):
    pdb=open(infile,'r')
    x=[0,]
    y=[0,]
    z=[0,]
    resnum=[0,]
    names=['empty',]
    for line in pdb:
        if not 'ATOM' and not 'HETATM 'in line:
            continue
        if (line[13:16]=='CA ') and (line[16]!='B') and (line[21]==chain):
            #print line
            x.append(float(line[28:38]))
            y.append(float(line[38:46]))
            z.append(float(line[46:54]))
            names.append(line[21:26])
            resnum.append(float(line[22:26]))
    x=np.asarray(x)
    y=np.asarray(y)
    z=np.asarray(z)
    resnum=np.asarray(resnum,dtype="int")
    maxres=np.max(resnum)
    #print 'The largest residue number I found was',maxres
    #print 'And there are',len(names),'residues in the structure...'
    C=np.empty((int(maxres+1),3))
    C[:]=np.NAN
    C[resnum,0]=x
    C[resnum,1]=y
    C[resnum,2]=z

    fullnames=[None]*(maxres+1)
    for n in range(len(resnum)):
        res=resnum[n]
        fullnames[res]=names[n]

    #print 'My fullnames matrix is size',len(fullnames)
    pdb.close()

    return C,fullnames

# Calculate a distance matrix
# This is the part where the distance between two 3D points are being calculated.
def distmat(C):
    x=np.power(diffmat(C[:,0]),2)
    y=np.power(diffmat(C[:,1]),2)
    z=np.power(diffmat(C[:,2]),2)
    D=np.power(x+y+z,0.5)
    return D

# calculate a difference matrix
# I think distance is being calculated and output as d
def diffmat(r):
    d=np.zeros((len(r),len(r)))
    for n in range(len(r)):
        for m in range(len(r)):
            d[n,m]=r[n]-r[m]
    return d

# extract unit cell from pdb file
# a function is defined, but I don't think it's called anywhere
def unitcell(pdbfile):
    infile=open(pdbfile,'r')
    for line in infile:
        elements=line.split()
        if elements[0]=='CRYST1':
            a=float(elements[1])
            b=float(elements[2])
            c=float(elements[3])
            alpha=float(elements[4])
            beta=float(elements[5])
            gamma=float(elements[6])
    return a,b,c,alpha,beta,gamma




pdb1,names1=readpdb(moleculeA,chainA)


#print 'I read',len(names1),'atoms..'




pdbs=[moleculeB]


labels=['']


for p in range(len(pdbs)):
    pdb2,names2=readpdb(pdbs[p],chainB)


    #print 'I read',len(names2),'atoms...'
    #print names1==names2
    label=labels[p]
    plotout='C:/Users/brain/Desktop/'+label+'.png'
    #print np.shape(pdb1),np.shape(pdb2)

    r=distmat(pdb2)-distmat(pdb1)
    a=distmat(pdb1)
    b=distmat(pdb2)

    r[0,:]=0
    r[:,0]=0
    mask=np.transpose(np.tri(r.shape[0],k=-1))
    r=np.ma.array(r,mask=mask)
    fig,ax=plt.subplots()

    cax=ax.matshow(r,cmap='bwr',clim=(-1.*maxdifference,maxdifference),interpolation='none')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


    l=[]
    s=[]
    e=[]
    sa=3
    st=5
    for n in range(len(l)):
        m=((s[n]+e[n])/2)-1.5*len(l[n])
        ax.text(m+st,m-st,l[n],rotation=-45)
        ax.annotate("",xy=(s[n]+sa,s[n]-sa),xycoords='data',xytext=(e[n]+sa,e[n]-sa),textcoords='data',arrowprops=dict(arrowstyle="<->",
                                connectionstyle="arc3"),)
    #ax.text(100,20,label, size=18)
    #ax.text(280,-10,r"$\Delta$$\vert$$\Delta$$\vec r$$\vert$ [$\AA$]",size=18)
    #ax.text(340,0,r"  $\blacktriangleleft$ $\blacktriangleright$")
    #ax.text(330,10,"  further apart",size=10)
    #ax.text(340,270,r"  $\blacktriangleright$$\blacktriangleleft$")
    #ax.text(330,280,"  closer together",size=10)

    cbar=fig.colorbar(cax)
    plt.savefig(plotout,dpi=300)
    plt.show()

