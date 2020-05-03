#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 00:58:24 2020

@author: jguterl
"""

class UEDGEMesh():
    def __init__(self):
        pass
    
    def PlotSeparatrix(self,ax=None,color='r',linewidth=1,**kwargs):
        from uedge import com
        sepx=np.concatenate((com.rm[:,com.iysptrx,3],np.array([com.rm[-1,com.iysptrx,4]])))
        sepy=np.concatenate((com.zm[:,com.iysptrx,3],np.array([com.zm[-1,com.iysptrx,4]])))
        
        if ax is None:
            ax=plt.gca()
        ax.plot(sepx,sepy,color=color,linewidth=linewidth,**kwargs)


    def ShowCell(self,ixiy,rm=None,zm=None,ax=None):
        if ax is None:
            ax=plt.gca()
        
        if type(ixiy)!=list:
            ixiy=[ixiy]
        if rm is None:
            rm=com.rm
        if zm is None:
            zm=com.zm
        for (ix,iy) in ixiy:
            r=rm[ix,iy,0]
            z=zm[ix,iy,0]
            r0 = [rm[ix, iy, 1], rm[ix, iy, 2],rm[ix, iy, 4], rm[ix, iy, 3], rm[ix, iy, 1]]
            z0 = [zm[ix, iy, 1], zm[ix, iy, 2], zm[ix, iy, 4], zm[ix, iy, 3], zm[ix, iy, 1]]
            ax.plot(r0, z0, 'r-', label='Grid', linewidth=2)
            
            annot = ax.annotate("ix={},iy={}".format(ix,iy), xy=(r,z), xytext=(-20,20),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))

    def PlotMesh(self,r,z,ax=None,Fill=False,Verbose=False,facecolor=None,edgecolor='black',Title=''):
        """Plot UEDGE grid."""
        if ax is None:
            ax=plt.gca()
            
        def onpick(evt):
            if evt.artist in Pos.keys():
                annot.set_visible(False)
                annot.xy = Pos[evt.artist]
                evt.artist.set_facecolor='blue'
                evt.artist.set_fill=True
                annot.set_text(Dic[evt.artist])
                annot.set_visible(True)
            if evt.mouseevent.button == 3:
                annot.set_visible(False)    
            ax.figure.canvas.draw_idle()
    
        Nx=len(r)
        Ny=len(r[0])
        ax.figure.suptitle("{} nx={} ny={}".format(Title,Nx,Ny))
        patches=[]
        if ax is None:
            ax=plt.gca()
        idx=[np.array([1,2,4,3,1])]
        Dic={}
        Pos={}
        Obj={}
        for i in range(Nx):
            for j in range(Ny):
                    Data=np.concatenate((r[i][j][idx],z[i][j][idx])).reshape(2,5).T
                    p=matplotlib.patches.Polygon(Data,closed=True,fill=False,edgecolor=edgecolor,label='ix={},iy={}'.format(i,j),picker=5)
                    if not Fill:
                        c=ax.add_patch(p) #draw the contours
                    else:
                            patches.append(p)
                            
                    Dic[p]='ix={},iy={}'.format(i,j)
                    Pos[p]=(r[i][j][0],z[i][j][0])
                    Obj[p]=c
    
        if Fill:
            Collec=matplotlib.collections.PatchCollection(patches,match_original=True,facecolors=None,edgecolors=edgecolor)
            c=ax.add_collection(Collec)
            
        plt.ylim(z.min(),z.max())
        plt.xlim(r.min(),r.max())
        annot = ax.annotate("", xy=(0,0), xytext=(-20,20),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        
        #Create event handler
        
            
        ax.figure.canvas.mpl_connect('pick_event', onpick)


    def ImportMesh(self,fname:str = 'gridue')->dict:
        """Import UEDGE grid file as dictionary."""
        try:
            f= open(fname, mode = 'r')
            Values = [int(x) for x in next(f).split()]
            HeaderItems = ['nxm', 'nym', 'ixpt1', 'ixpt2', 'iyseptrx1']
            gridue_params=dict(zip(HeaderItems, Values))
            next(f)
            BodyItems = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
            Str={ i : [] for i in BodyItems }
            k=iter(Str.keys())
            Key=next(k)
            for line in f:
                if line=='iogridue\n':
                    continue
                if line=='\n':
                    try:
                        Key=next(k)
                    except:
                        continue
                    print(Key)
                else:
                    Str[Key].append(line)
            f.close()
            nx=gridue_params['nxm']+2
            ny=gridue_params['nym']+2
            for k,v in Str.items():
                L=(''.join(v).replace('\n','').replace('D','e')).split()
                l=iter(L)
                vv= next(l)
    
                data_=np.zeros((nx,ny,5))
                for n in range(5):
                        for j in range(ny):
                            for i in range(nx):
    
                                data_[i][j][n]=float(vv)
    
                                try:
                                    vv=next(l)
                                except:
                                    continue
                gridue_params[k]=data_
            return gridue_params
        except Exception as e:
            print(repr(e))