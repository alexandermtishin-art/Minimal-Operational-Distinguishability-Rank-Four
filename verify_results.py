"""Exact checks for the explicit examples; proofs remain in the manuscripts.
Run: python3 verify_results.py
No third-party dependencies and no floating-point arithmetic.
"""
from fractions import Fraction as F
from itertools import product, permutations
import json

def mul(a,b):
 return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))) for i in range(len(a)))
def check():
 results={}
 z=((0,0,0),(0,0,0),(0,0,1))
 n=((F(1,2),0,F(1,2)),(0,0,0),(F(1,2),0,F(1,2)))
 assert mul(z,z)==z and mul(n,n)==n and mul(z,n)!=mul(n,z)
 assert tuple(row[2] for row in mul(z,n))==(0,0,F(1,2))
 assert tuple(row[2] for row in mul(n,z))==(F(1,2),0,F(1,2))
 results['dephasing']='Exact rational Bloch matrices reproduce Example 3.3'
 P=((1,0),(0,0));Q=((F(1,2),F(1,2)),(F(1,2),F(1,2)))
 assert mul(P,P)==P and mul(Q,Q)==Q and mul(P,Q)!=mul(Q,P)
 results['filters']='Idempotence and noncommutation verified'
 S={(0,0):set(product(range(2),repeat=2)),(0,1):{(0,0),(1,0),(1,1)},(1,0):{(0,0),(0,1),(1,1)},(1,1):{(0,1),(1,0),(1,1)}}
 G=[v for v in product(range(2),repeat=4) if all((v[x],v[2+y]) in S[x,y] for x,y in S)]
 assert len(G)==5 and not any(v[0]==v[2]==0 for v in G)
 assert all({a for a,b in s}=={0,1} and {b for a,b in s}=={0,1} for s in S.values())
 results['hardy_support']='Five global assignments; none extends 00 at context 00; all local projections full'
 # Reset witness: normal-form composition evaluated independently as maps.
 def act(f,s):
  m,A=f;k,B=s
  return k+m,(B|A if m==0 else A)
 def comp(f,g):
  m,A=f;n,B=g
  return (n,A|B) if m==0 else (m+n,A)
 forms=list(product(range(3),range(8)))
 for f,g,s in product(forms,forms,product(range(2),range(8))):assert act(comp(f,g),s)==act(f,act(g,s))
 assert [A for A in range(8) if sum(A&(1<<i)!=0 for i in range(3)) in (0,3)]==[0,7]
 results['reset_witness']='Normal-form multiplication and the two frame-neutral Boolean states verified'
 # All nonempty relations on a 3-by-2 set; all bijective local relabellings.
 cells=list(product(range(3),range(2)))
 def rectangular(D):return len(D)==len({a for a,b in D})*len({b for a,b in D})
 count=0
 for mask in range(1,1<<len(cells)):
  D={p for i,p in enumerate(cells) if mask>>i&1}
  for f,g in product(permutations(range(3)),permutations(range(2))):
   E={(f[a],g[b]) for a,b in D}
   assert rectangular(D)==rectangular(E);count+=1
 results['local_injective_preservation']={'finite_cases':count}
 # Free-product normal forms: atoms i=1,2,4; t represented by ('t',1).
 def norm(word):
  out=[]
  for typ,val in word:
   if val==0:continue
   if out and out[-1][0]==typ:
    prev=out.pop()[1];out.append((typ,prev|val if typ=='b' else prev+val))
   else:out.append((typ,val))
  return tuple(out)
 letters=[(('b',1),),(('b',2),),(('b',4),),(('t',1),)]
 words={()}
 for length in range(1,6):
  for inds in product(range(4),repeat=length):words.add(norm(sum((letters[i] for i in inds),())))
 def image(w):
  r=((F(1),F(0),F(0)),(F(0),F(1),F(0)),(F(0),F(0),F(1)))
  for typ,val in w:r=mul(r,z if typ=='b' else n)
  return r
 for w in words:
  for l in letters:
   assert image(norm(w+l))==mul(image(w),image(l))
 assert all(image((('t',m),))==n for m in range(1,20))
 results['free_product']={'normal_forms_through_length_5':len(words),'checked':'Generator relations, homomorphism onto the dephasing monoid of Example 3.3 (x_i -> D_z, t -> D_n) and distinct-degree powers'}
 return results
if __name__=='__main__': print(json.dumps(check(),indent=2))
