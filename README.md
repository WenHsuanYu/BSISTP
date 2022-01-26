# BSISTP
Balanced Selected-Internal Steiner Tree

I use approximate algorithms to find approximate solutions to BSISTP.

## Conclusion Algorithm 4 Time Complexity: (O(E^2)):

1. Use algorithm 1 for the BoSTP(bottleneck Steiner tree problem) to find a Steiner tree S for R in G.
2. Use algorithm 2 to find a balanced Steiner tree S for R in G.
3. If S is not a Selected-internal Steiner tree for R and R' in G 
   then
	3.1 Select any two vertices{𝑣𝑎, 𝑣𝑏} in R \ R'.
	3.2 Apply algorithm 3 to find a Selected-internal Steiner tree T in S^３.
---

## Algorithm 1 Time Complexity: O(E log E):
KRUSKAL(G):<br>
A = ∅<br>
For each vertex v ∈ G.V:<br>
	&nbsp;MAKE-SET(v)<br>
For each edge (u, v) ∈ G.E ordered by increasing order by weight(u, v):<br>
	&nbsp;if FIND-SET(u) ≠ FIND-SET(v):<br>
	&nbsp;A = A ∪ {(u, v)}<br>
	&nbsp;UNION(u, v)<br>
return A<br>

---
## Algorithm 2 Time Complexity: O(|E|^2):

![](https://github.com/WenHsuanYu/BSISTP/blob/main/pic/alg2.png)

---
## Algorithm 3 Time Complexity: O(|V|^2)
![](https://github.com/WenHsuanYu/BSISTP/blob/main/pic/alg3.png)

