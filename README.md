# BSISTP
Balanced Selected-Internal Steiner Tree

I use approximate algorithms to find approximate solutions to BSISTP.

## Conclusion Algorithm 4 (**O(E^2)**):

1. Use algorithm 1 for the BoSTP(bottleneck Steiner tree problem) to find a Steiner tree S for R in G.
2. Use algorithm 2 to find a balanced Steiner tree S for R in G.
3. If S is not a Selected-internal Steiner tree for R and R' in G 
   then
	3.1 Select any two vertices{𝑣𝑎, 𝑣𝑏} in R \ R'.
	3.2 Apply algorithm 3 to find a Selected-internal Steiner tree T in S^３.
---

## Algorithm 1
KRUSKAL(G):
A = ∅
For each vertex v ∈ G.V:
    MAKE-SET(v)
For each edge (u, v) ∈ G.E ordered by increasing order by weight(u, v):
    if FIND-SET(u) ≠ FIND-SET(v):       
    A = A ∪ {(u, v)}
    UNION(u, v)
return A

---
## Algorithm 2 Time Complexity: O(|E|^2)

![](https://github.com/WenHsuanYu/BSISTP/blob/main/pic/alg2.png)

---
## Algorithm 3
![](https://github.com/WenHsuanYu/BSISTP/blob/main/pic/alg3.png)

