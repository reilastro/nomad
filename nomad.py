#!/usr/bin/env python

import os
import sys

import math
R2D=180.0/math.pi
D2R=math.pi/180.0

from pprint import pprint

import time
import numpy as np
import scipy
from scipy import weave
from scipy.weave import converters

import random

try:
    import matplotlib.pyplot as plt
    NOMPL=False
except:
    NOMPL=True

def benchmark():
    N = 5000

    A = np.random.rand(N, N)
    B = np.random.rand(N, N)
    C = np.zeros([N, N], dtype=float)

    t = time.clock()
    weave_inline_loop(A, B, C, N)
    print time.clock() - t


def weave_inline_loop(A, B, C, N):
    code = """
           int i, j;
           for (i = 0; i < N; ++i)
           {
               for (j = 0; j < N; ++j)
               {
                   C(i, j) = A(i, j) * B(i, j);
               }
           }
           """
    weave.inline(code, ['A', 'B', 'C', 'N'], type_converters=converters.blitz, compiler='gcc')

class nomad():
    def __init__(self,location):
        self.all_cat=None     # A list of all catalog files
        self.record_size=22   # 22 integers in catalog file
        self.byte_size=4      # "i4" specifically
        self.dir=location
        self.decbin=0.1
        self.rabin=15.0
        self.data=self.no_stars()
        # Integers in cat files to units
        """
   ( 1)   RA at 2000.0 in integer 0.001 arcsec
   ( 2)   SPD at 2000.0 in integer 0.001 arcsec
   ( 3)   std. dev. of RA*COS(dec) in integer 0.001 arcsec at central epoch
   ( 4)   std. dev. of SPD in integer 0.001 arcsec at central epoch
   ( 5)   proper motion of RA*COS(dec) in integer 0.0001 arcsec/year
   ( 6)   proper motion of SPD in integer 0.0001 arcsec/year
   ( 7)   std. dev. of (5) in integer 0.0001 arcsec/year
   ( 8)   std. dev. of (6) in integer 0.0001 arcsec/year
   ( 9)   central epoch of RA in integer 0.001 year
   (10)   central epoch of SPD in integer 0.001 year
   (11)   B magnitude in integer 0.001 mag
   (12)   V magnitude in integer 0.001 mag
   (13)   R magnitude in integer 0.001 mag
   (14)   J magnitude in integer 0.001 mag
   (15)   H magnitude in integer 0.001 mag
   (16)   K magnitude in integer 0.001 mag
   (17)   USNO-B1.0 ID integer
   (18)   2MASS ID integer
   (19)   YB6 ID integer
   (20)   UCAC-2 ID integer
   (21)   Tycho2 ID integer
   (22)   flags integer
"""

        mult=scipy.ones((self.record_size),dtype=float)
        try:
            mult[0:4]*=0.001 # RA, DEC, sig RA, sig DEC
            mult[0:2]/=3600.0 # degrees
            mult[4:8]*=0.0001 # proper motions and sigmas
            mult[8:10]*=0.001 # epoch os RA DEC
            mult[10:16]*=0.001 # Mags BVRJHK
            self.load_mult=mult
        except:
            print "Unable to set load_mult[iplier]"            
        self.load_units=['deg','deg',
                         'arcsec','arcsec',
                         'arcsec/yr','arcsec/yr',
                         'arcsec/yr','arcsec/yr',
                         'yr','yr',
                         'mag','mag','mag','mag','mag','mag',
                         'int','int','int','int','int',
                         'int flag',
                         ]
        self.load_vars=['RA','DEC','sRA','sDEC',
                        'RAmotion','DECmotion','sRAmotion','sDECmotion',
                        'RAepoch','DECepoch',
                        'B','V','R','J','H','K',
                        'flag']
        if not os.path.exists(self.dir):
            raise Exception("Nomad location does not exist")
        pass

    def get_dir(self):
        return self.dir
    def set_dir(self,location):
        if os.path.exists(location):
            self.dir=location
        else:
            print "Bad nomad location",location
        return self.dir

    def get_bl(self,ra=None,dec=None):
        """
        http://scienceworld.wolfram.com/astronomy/GalacticCoordinates.html
        """
        if ra==None: ra=self.get_column("RA")
        if dec==None: dec=self.get_column("DEC")
        if type(ra)==float: 
            ral=scipy.zeros(1)
            ral[0]=ra
            ra=ral
        if type(dec)==float: 
            decl=scipy.zeros(1)
            decl[0]=dec
            dec=decl

        c62=math.cos(62.6*D2R)
        s62=math.sin(62.6*D2R)
        
        b=scipy.sin(dec*D2R)*c62
        b-=scipy.cos(dec*D2R)*scipy.sin((ra-282.25)*D2R)*s62
        b=scipy.arcsin(b)*R2D
        
        cosb=scipy.cos(b*D2R)
        #l minus 33 degrees
        lm33=(scipy.cos(dec*D2R)/cosb)*scipy.cos((ra-282.25))
        l=scipy.arccos(lm33)*R2D+33.0
        return b,l
        



    def load_acc(self,fname):
        # File is
        # Start hour, star1 (1 base), Nstar
        # Rval is
        # start degree, star1 (1 base), Nstar

        rval=[]
        i=0
        for line in open(fname):
            i+=1
            l=line.strip()
            l2=l.split()
            rval.append([float(l2[0])*15.0,int(l2[1]),int(l2[2])])
        return rval

    def all_cat_files(self):
        if self.all_cat==None:
            self.all_cat=[]
            for i in range(-900,900):
                i2="%04d"%(900+i)
                f=os.path.join(self.dir,i2[:3],'m'+i2+'.cat')
                self.all_cat.append(f)
        return self.all_cat

    def angular_distance(self,ra1,dec1,ra2,dec2):
        cosa=math.sin((90.0-dec1)*D2R)*math.sin((90.0-dec2)*D2R)
        cosa*=math.cos((ra1-ra2)*D2R)
        cosa+=(math.cos((90.0-dec1)*D2R)*math.cos((90.0-dec2)*D2R))
        return math.acos(cosa)*R2D

    def load_cat(self,fname,maxsize=1E5,**kwargs):
        current_dec=self.file_to_dec(fname)
        mindec=-90
        maxdec=90
        minra=0
        maxra=360
        
        # Skip dec's
        if 'ra' in kwargs and 'dec' in kwargs and 'radius' in kwargs:
            ra=kwargs['ra']
            dec=kwargs['dec']
            radius=kwargs['radius']
            mindec=dec-radius
            maxdec=dec+radius

        if 'mindec' and 'maxdec' in kwargs:
            mindec=kwargs['mindec']
            maxdec=kwargs['maxdec']

        if current_dec<mindec-self.decbin*0.5: 
            #print "{", ; sys.stdout.flush()
            #print "Skip",current_dec
            return self.no_stars()
        if current_dec>=maxdec+self.decbin*0.5:
            #print "}", ; sys.stdout.flush()
            #print "Skip",current_dec
            return self.no_stars()



        if 'ra' in kwargs and 'dec' in kwargs and 'radius' in kwargs:
            ra=kwargs['ra']
            dec=kwargs['dec']
            radius=kwargs['radius']
            mindec=dec-radius
            maxdec=dec+radius
            try:
                deltara1=math.acos(math.fabs(radius*math.cos(D2R*current_dec-self.decbin*0.5)))
            except:
                deltara1=0.0
            try:
                deltara2=math.acos(math.fabs(radius*math.cos(D2R*current_dec+self.decbin*0.5)))
            except:
                deltara2=0.0
            deltara=max((deltara1,deltara2))
            minra=ra-deltara
            maxra=ra+deltara+self.rabin
            if minra<0 and maxra>360:
                minra=0.0
                maxra=360
            else:
                if maxra>360: maxra-=360
                if minra<0: minra+=360

        if 'minra' and 'maxra' in kwargs:
            minra=kwargs['minra']
            maxra=kwargs['maxra']


        rval=self.no_stars()
        try:
            accel=self.load_acc(fname[:-4]+".acc")
        except:
            print "No accel file for",fname
            return self.no_stars()
        i=0
        f=open(fname,'rb')

        for chunk in accel:
            startra=chunk[0]
            if startra+self.rabin<minra: 
                #print "-", ; sys.stdout.flush()
                #print "Skip",chunk
                continue
            if startra>maxra: 
                #print "-", ; sys.stdout.flush()
                #print "Skip",chunk
                continue
            #print "+", ; sys.stdout.flush()
        
            start=self.record_size*(chunk[1]-1)*self.byte_size
            f.seek(start,0)
            nstar=chunk[2]


            chunkload=scipy.fromfile(f,dtype='i4',count=self.record_size*nstar)
            chunkload=scipy.reshape(chunkload,(nstar,self.record_size))*self.load_mult
            chunkload[:,1]-=90.0
            #*self.load_mult
            #print scipy.shape(chunkload)
            #print scipy.shape(chunkload[:,0])
            rval=scipy.append(rval,chunkload,axis=0)
            #print i,chunk[0],chunk,len(chunkload),len(rval)
            #print chunk,scipy.shape(rval)
        #print
        return rval

    def no_stars(self):
        return scipy.empty((0,22))

    def get_cut_radius(self,ra,dec,radius,declist=None,ralist=None,data=None):
        # Cut values
        cra=ra
        cdec=dec
        cradius=radius
        # Arrays from star list
        if declist==None:
            declist=self.get_column('DEC')
        if ralist==None:
            ralist=self.get_column('RA')

        if declist==None or len(declist)<1:
            return None
        if ralist==None or len(ralist)<1:
            return None

        radius=self.get_radius(ra=cra,dec=cdec,
                               ralist=ralist,declist=declist)
        cut=radius<cradius
        return cut

    def get_radius(self,ra,dec,ralist=None,declist=None):
        """ ra is single value of RA in degrees from which radius is measured
        dec is is single value of DEC in degrees from which radius is measured
        ralist is the 1d array of RA for stars to consider 
        declist is the 1d array of DEC for stars to consider 
        if a list of ra or dec are None then extract from self.data using get_column
        """
        # Arrays from star list
        if declist==None:
            declist=self.get_column('DEC')
        if ralist==None:
            ralist=self.get_column('RA')

        # Anuglar distance from ra,dec for all stars
        cosa=scipy.sin((90.0-declist)*D2R)*scipy.sin((90.0-dec)*D2R)
        cosa*=scipy.cos((ralist-ra)*D2R)
        cosa+=scipy.cos((90.0-declist)*D2R)*scipy.cos((90.0-dec)*D2R)
        # radius array
        radius=scipy.arccos(cosa)*R2D
        return radius

    def get_column(self,k):
        try:
            i=self.load_vars.index(k)
        except:
            return None
        try:
            return self.data[:,i]
        except:
            return None

    def load_cats(self,**kwargs):
        rval=self.no_stars()

        for f in self.all_cat_files():
            new_stars=self.load_cat(f,**kwargs)
            if new_stars==None or len(new_stars)==0: continue
            if rval==None: rval=new_stars
            else: rval=scipy.append(rval,new_stars,axis=0)

        self.data=rval
        #print "R"
        if 'ra' in kwargs and 'dec' in kwargs and 'radius' in kwargs:
            cut=self.get_cut_radius(ra=kwargs['ra'],
                                    dec=kwargs['dec'],
                                    radius=kwargs['radius'])
            data=self.apply_cut(cut)
            self.set_data(data)
        return rval




    def dec_to_file(self,dec):
        # - NOMAD id numbers consist of a zone number (0 to 1799
        # in 1/10 degree steps of declination, starting from the
        # south celestial pole), and a running star number in that
        # zone (increasing along RA).
        dec2="m%04d"%(int(math.floor((dec+90)*10)))
        f=os.path.join(dec2[1:4],dec2)
        rval={}
        rval['acc']=os.path.join(self.dir,f+".acc")
        rval['cat']=os.path.join(self.dir,f+".cat")
        for k in rval:
            if not os.path.exists(rval[k]):
                print "No",k,rval[k]
        print dec,rval['cat'],self.file_to_dec(rval['cat'])
        return rval

    def file_to_dec(self,f):
        # - NOMAD id numbers consist of a zone number (0 to 1799
        # in 1/10 degree steps of declination, starting from the
        # south celestial pole), and a running star number in that
        # zone (increasing along RA).
        return int(os.path.basename(f).split(".")[0][-4:])/10.0-90.0+self.decbin*0.5

    def get_data(self):
        return self.data
    def set_data(self,data):
        self.data=data

    def get_cut(self, k, val, oper='>',col=None):
        """ 
        cut on column k using operator oper
        if col is None use self.get_column k
        """
        if col==None:
            try:
                col=self.get_column(k)
            except:
                print "Cannot cut on var",k
                return
        if oper=='>': cut=col>val
        if oper=='<': cut=col<val
        if oper=='=': cut=col==val
        return cut

    def apply_cut(self,cut,data=None):
        """ Cut is an array of True False of length data
        if data==None is self.data
        """
        #self.data=scipy.delete(self.data,cut)
        if data==None:
            data=self.get_data()
            if len(data)==0: return self.no_stars()
            data=data[cut]
            return data
        else:
            if len(self.data)==0: return self.no_stars()
            self.data=self.data[cut]
            return self.data

    def plot(self,**kwargs):
        if self.data==None: 
            print "No data loaded"
            return
        #print kwargs
        data={}
        for k in ['x','y','z','s','c']:
            if k in kwargs:
                try:
                    if type(kwargs[k])==type('string'):
                        data[k]=self.data[:,self.load_vars.index(kwargs[k])]
                    else:
                        data[k]=kwargs[k]
                    print k,len(data[k]),data[k]
                except:
                    print "Error in plot var",k

        #print data.keys(),data['x'][:4]
        if 'x' in data.keys() and 'y' in data.keys():
            plt.scatter(**data)
        elif 'x' in data:
            data['bins']=50
            plt.hist(**data)
        plt.show()

    def load_raw(self,ifile,nstar=0):
        f=open(ifile,"rb")
        data=f.read(self.record_size*self.byte_size*nstar)
        f.close()
        return data
        
    def rewrite_file(self,ifile,odir,cut=None):
        ofile=ifile.replace(self.dir,odir)
        #if os.path.exists(ofile):
        #    print "%s exists"%ofile
        #    return
        odir=os.path.dirname(ofile)
        if not os.path.exists(odir):
            try:
                os.mkdir(odir)
            except:
                print "Failed to create",odir
        try:
            accel=self.load_acc(ifile[:-4]+".acc")
        except:
            print "No accel file for",ifile
            return 
        nstar=accel[-1][1]+accel[-1][2]-1
        data=self.load_raw(ifile,nstar)
        ra=self.get_column('RA')
        R=self.get_column('R')


        odata=""
        for i in range(nstar):
            if not cut[i]: continue
            odata+=data[i*self.record_size*self.byte_size:i*self.record_size*self.byte_size+self.record_size*self.byte_size]

        cbin_start_rah=0.00
        nbin_start_rah=0.25
        cbin_start_star=1
        cbin_n_star=0
        chunks=[]
        for i in range(nstar):
            if not cut[i]: continue
            while ra[i]>nbin_start_rah*15 and cbin_start_rah<24.00:
                chunks.append([cbin_start_rah,cbin_start_star,cbin_n_star])
                
                cbin_start_rah=nbin_start_rah
                nbin_start_rah=cbin_start_rah+0.25
                cbin_start_star+=cbin_n_star
                cbin_n_star=0
            cbin_n_star+=1
        while (len(chunks)<24*4):
            chunks.append([cbin_start_rah,cbin_start_star,cbin_n_star])
                
            cbin_start_rah=nbin_start_rah
            nbin_start_rah=cbin_start_rah+0.25
            cbin_start_star+=cbin_n_star
            cbin_n_star=0


        fout=open(ofile,"wb")
        fout.write(odata)
        fout.close()
        print "   Wrote",
        print "%d"%(len(odata)),
        print "%d"%(len(odata)/self.byte_size),
        print "%d"%(len(odata)/self.byte_size/self.record_size),
        print

        fout=open(ofile[:-4]+".acc","w")
        for chunk in chunks:
            oline="%5.2f %11d %11d\n"%(chunk[0],chunk[1],chunk[2])
            fout.write(oline)
            #print oline,
        fout.close()
        #print "   Wrote %r"%(chunks[0]),
        #print "   Wrote %r"%(chunks[-1])



        pass

def main():
    if sys.platform=='darwin' or os.path.exists("/nomad"):
        s=nomad("/nomad")
    else:
        s=nomad("/nfs/slac/g/ki/ki19/des/nomad")

    print "Nomad loading from",s.get_dir()


    if 1==2:
        flist=s.all_cat_files()
        flist.reverse()
        for f in flist:
            print f
            data=s.load_cat(f)
            print "   Loaded %d stars"%(len(data))
            s.set_data(data)
            b,l=s.get_bl()
            cut1=scipy.fabs(b)>30
            cut2=s.get_cut('R',15,'>')
            cut3=s.get_cut('R',19,'<')
            cut=cut1*cut2*cut3
            print "   Writing %d stars"%((cut1*cut2*cut3).sum())
            s.rewrite_file(f,odir="/nomad2",cut=cut1*cut2*cut3)

        return


    
    if 1==1:
        print "Nomad loading from",s.get_dir()
        N=200
        n=scipy.zeros(N)
        blist=scipy.zeros(N)
        i=-1
        while (i<N-1):
            i+=1
            dec=random.randrange(-20,90)
            ra=random.randrange(0,360)
            b,l=s.get_bl(ra=ra,dec=dec)
            #print "DEC=%6.1f RA=%7.1f     b=%5.1f     l=%5.1f"%(dec,ra,b,l)
            if math.fabs(b)<30: 
                i-=1
                continue
            
            tot=0
            curdec=dec
            curra=ra
            # Load list of stars at RA, DEC
            s.load_cats(dec=curdec,ra=curra,radius=0.107)
            #print "#1",i,len(s.get_data())
            
            # Cut on R>15
            cut=s.get_cut('R',15,'>')
            data=s.apply_cut(cut)
            s.set_data(data)
            #print "#2",i,len(s.get_data())

            # Cut on R<19
            cut=s.get_cut('R',19,'<')
            data=s.apply_cut(cut)
            s.set_data(data)
            #print "   #3",i,len(s.get_data())

            tot+=len(s.get_data())
            if 1==2:
                # Array of circle sizes to draw stars
                star_size=5.0+scipy.power(19-s.get_column('R'),3)
                # Draw stars
                s.plot(x="RA",y='DEC',s=star_size)
            n[i]=tot
            blist[i]=b
            print "Total",i,b,tot
        s.plot(y=n,x=blist)
        return s
    return

    for i in range(10):
        # Active area of 1 GFA is 0.006 deg sq
        # Active area of 6 GFA is 0.036

        # Total area is equal to a circle with radius 0.107 degrees
        # Pi r^2 = 3 * 0.1 * 0.1 = 0.03

        ra=random.randrange(0,360)
        dec=random.randrange(45,90)

        # Load list of stars at RA, DEC
        s.load_cats(dec=dec,ra=ra,radius=0.107)

        # Array of circle sizes to draw stars
        star_size=5.0+scipy.power(19-s.get_column('R'),3)
        # Draw stars
        s.plot(x="RA",y='DEC',s=star_size)
        break
    return s
    for dec in range(-90,91,10):
        files=s.dec_to_file(dec)
        #s.load_acc(files['acc'])
        s.load_cat(files['cat'])
        print dec,s.dec_to_file(dec)
        #break

if __name__=="__main__":
    s=main()
    #benchmark()
