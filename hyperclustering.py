import re
#import scipy
from pylab import *
import numpy as np
from numpy import *
import random

def random_assignment(K,M):
    return [random.randint(0,K-1) for i in range(0,M)]


def distance(point,hyp):
	#assert( len(point) == (len(hyp)-1) > 0 )
	if not (len(point) == (len(hyp)-1) > 0):
		return -1
	dist=abs( dot(hyp[0:-1],point) - hyp[-1])
	dist=dist/linalg.norm(hyp[0:-1]) #sqrt(dot(hyp[0:-1],hyp[0:-1]))
	return dist

def reassign_points(points, hyperplanes):
    
    tot_dist=0
    assignment=[0 for x in range(len(points))]
    
    for j in range(len(points)):
        min_dist= +infty
        min_tag=0
        for i in range(len(hyperplanes)):
            #print points[j],
            #print hyperplanes[i]
            d=distance(points[j],hyperplanes[i])
            #TODO can be done better?
            if d<min_dist:
                min_dist=d
                min_tag=i
        assignment[j]=min_tag
        tot_dist+=min_dist
    return assignment
                    
def tot_distance(points, hyperplanes,assign):
    tot_dist=0
    for j in range(len(points)):
         tot_dist+=distance(points[j],hyperplanes[assign[j]])
    return tot_dist
                         


def hclust(points,N,K=3,max_it=1000):
    #K idem
  
	#assert(reduce( (lambda x,y : x if (x==y) else -1), map(len, points) ) == N)
    
	M = len(points) #numero di punti
	assignment = random_assignment(K,M)
	print assignment
	halting=False
	iteration=0
	hyperplanes=[[] for x in range(K)]
	old_dist=f_min=+infty
	best_assign=[]
	
	while iteration<max_it:
		for j in range(K):
			Aj = matrix(compress(equal(assignment, j), points, 0))
			hyperplanes[j]=update(Aj)

		
		new_dist= tot_distance(points, hyperplanes, assignment)
		
		
		print 'F:',
		print new_dist
        
		#if tot_dist >= f_min:
		#	break
		
		if assignment==best_assign or new_dist==old_dist:
			break
		
		if new_dist < f_min:
			f_min=new_dist
			best_assign=assignment
			best_hyp=hyperplanes
		
		old_dist=new_dist

		assignment=reassign_points(points, hyperplanes)
		iteration+=1
	
	print iteration, new_dist
	return best_assign, best_hyp
	

def update(Aj):
	if len(Aj)==0:
		return [random.random() for i in range(Aj.ndim+1)]
	if len(Aj)==1:
		gamma=0
		hyperplane=[0 for x in range(Aj.ndim+1)]
		hyperplane[1]=1
		hyperplane[Aj.ndim]=Aj[0,1]
		return hyperplane
	
	B=Aj.T*((eye(len(Aj))-ones([len(Aj),len(Aj)])/len(Aj)))*Aj
	
	eigvals,eigvecs=linalg.eig(B)
	# w = min eig of A^t*B*A where B = [I - pp^T/p^Tp ], i.e., with p=[1]^m, B = [I_m - ones(m)/m]
	# gamma = p^t*A*w/p^T*p
	
	min_val=(eigvals[0])
	w=array(eigvecs.T[0])
	for i in range(1,len(eigvals)):
		if eigvals[i]<min_val:
			min_val=eigvals[i]
			w=array(eigvecs.T[i])
	e=matrix(ones([len(Aj),1]))
	gamma=(e.T*Aj*w.T/len(Aj)).tolist()[0][0]
	
	hyperplane=w.tolist()[0]
	hyperplane.append(gamma)    
	return hyperplane
	#the column eigvecs[:,i] is the eigenvector corresponding to the eigenvalue w[i]
	


def hplot(points,hyperplanes=[], assign=[]): #gamma sta alla fine. hyperplanes contengono w di dim N e gamma, totale=N+1
    if assign==[] and hyperplanes==[]:
        x = array([u[0] for u in points])
        y = array([u[1] for u in points])
        #xfit=array([amin(x), amax(x)])
        xfit=array([amin(x) - 0.1*(amax(x)- amin(x)), amax(x)+ 0.1*(amax(x)-amin(x))])
        #ion()
        plot(x,y,'ks') 
    elif assign==[]:
        x = array([u[0] for u in points])
        y = array([u[1] for u in points])
        #xfit=array([amin(x), amax(x)])
        xfit=array([amin(x) - 0.1*(amax(x)-amin(x)), amax(x)+ 0.1*(amax(x)-amin(x))])
        #ion()
        plot(x,y,'ks') 
        for k in range(0,len(hyperplanes)):
            fit=array((-hyperplanes[k][0]*xfit + hyperplanes[k][len(hyperplanes[k])-1])/hyperplanes[k][1])#a*x+b
            plot(xfit,fit,'r-', lw=2)

        show() #CASO CON = 0 (verticale)
    else:
        markers=['s','+' ,'s', '*' , ',' , '.' , '1' , '2' , '3' , '4' , '<' , '>' , 'D' , 'H' , '^' , '_' , 'd' , 'h' , 'o' , 'p' , 's' , 'v' , 'x' , '|']
        colors=['b','g','r','c','m','y','k','w']
        for j in range(len(hyperplanes)):
            p=compress(equal(assign, j), points, 0)
            if len(p)>=1:                
                x = array([u[0] for u in p])
                y = array([u[1] for u in p])
                xfit= (array([amin(x) - 0.1*(amax(x)-amin(x)), amax(x)+ 0.1*(amax(x)-amin(x))]) if amax(x)>amin(x) else array([x-0.1*x, x+0.1*x]))
                
                plot(x,y,colors[j%len(colors)]+markers[j%len(markers)])
                fit=array((-hyperplanes[j][0]*xfit + hyperplanes[j][len(hyperplanes[j])-1])/hyperplanes[j][1])#a*x+b
                plot(xfit,fit,colors[j%len(colors)]+'-', lw=2)
                for pi in p:
					if hyperplanes[j][0]!=0 and hyperplanes[j][1]!=0:
						ics=(hyperplanes[j][len(hyperplanes[j])-1]/hyperplanes[j][1]-pi[1]+hyperplanes[j][1]/hyperplanes[j][0]*pi[0])*hyperplanes[j][1]*hyperplanes[j][0]/(hyperplanes[j][1]*hyperplanes[j][1]+hyperplanes[j][0]*hyperplanes[j][0])
						xfit=array([min([ics,pi[0]]),max([ics,pi[0]])])
						fit=array(pi[1] + hyperplanes[j][1]/hyperplanes[j][0]*xfit - hyperplanes[j][1]/hyperplanes[j][0]*pi[0])#
						plot(xfit,fit,colors[j%len(colors)]+'--', lw=1)
            else:
                x = array([u[0] for u in points])
                y = array([u[1] for u in points])
                for pi in p:
                    plot(pi[0],pi[1],colors[j%len(colors)]+markers[j%len(markers)])
                xfit=array([amin(x) - 0.1*(amax(x)-amin(x)), amax(x)+ 0.1*(amax(x)-amin(x))])
                fit=array((-hyperplanes[j][0]*xfit + hyperplanes[j][len(hyperplanes[j])-1])/hyperplanes[j][1])#a*x+b
                plot(xfit,fit,colors[j%len(colors)]+'-', lw=2)

            
        show()
            
                
      
    

def load_points(file_name):
    f=open(file_name)
    points=[]
    for line in f:
        point=[]
        for num in line.split():
            try:
		point.append(float(num))
	    except:
                print 'not a float, skipped: '+num
		pass
	point= array(point)   
        points.append(point)
    return points

def main():
    
    points=load_points("/home/leo/prova")
    
    assignment, hyps=hclust(points,2,3)
    print assignment
    print hyps
    
    hplot(points,hyps, assignment)
    
    return 0

if __name__ == '__main__':
	main()

