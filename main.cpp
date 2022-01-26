#include <bits/stdc++.h>

using namespace std;

typedef pair<int,int> vertices;
typedef vector <pair<int,vertices>> weight_edge;
static const int vertex = 6;
static const int edges = 15;

struct Graph{
    int V,E;
    vector<pair<int,vertices>> edges;   

	Graph(int V, int E)
	{
		this->V = V;
		this->E = E;
	}   

    void addEdge(int u, int v, int w)
	{
		edges.push_back({ w, {u, v} });
	}

    vector <pair<int,vertices>> kruskalMST();
};

struct DisjointSets
{
	int *parent, *rnk;
	int n;

	// Constructor. 
	DisjointSets(int n)
	{
		// Allocate memory 
		this->n = n;
		parent = new int[n+1];
		rnk = new int[n+1];

		// Initially, all vertices are in 
		// different sets and have rank 0. 
		for (int i = 0; i <= n; i++)
		{
			rnk[i] = 0;

			//every element is parent of itself 
			parent[i] = i;
		}
	}

	// Find the parent of a node 'u' 
	// Path Compression 
	int find(int u)
	{
		/* Make the parent of the nodes in the path
		from u--> parent[u] point to parent[u] */
		
		if (u != parent[u])
			parent[u] = find(parent[u]);
		return parent[u];
	}

	// Union by rank 
	void merge(int x, int y)
	{
		x = find(x), y = find(y);

		/* Make tree with smaller height
		a subtree of the other tree */
		if (rnk[x] > rnk[y])
			parent[y] = x;
		else // If rnk[x] <= rnk[y] 
			parent[x] = y;

		if (rnk[x] == rnk[y])
			rnk[y]++;
	}
};

vector <pair<int,vertices>> Graph::kruskalMST()
{
    int mst_wt = 0; // Initialize result

    // Sort edges in increasing order on basis of cost
    sort(edges.begin(), edges.end());
  
    // Create disjoint sets
    DisjointSets ds(V);
	vector<pair<int,vertices>> tree; 
    // Iterate through all sorted edges
    vector< pair<int, vertices> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
  
        int set_u = ds.find(u);
        int set_v = ds.find(v);
  
        // Check if the selected edge is creating
        // a cycle or not (Cycle is created if u
        // and v belong to same set)
        if (set_u != set_v)
        {
  
            // Update MST weight
            mst_wt += it->first;
  
            // Merge two sets
            ds.merge(set_u, set_v);
			tree.push_back({ it->first, {u, v} });
        }
    }
	cout << "size of the tree: " << tree.size()<< endl;
    return tree;
}
//Time Complixity: O(|V||E|)
weight_edge sort_edge(Graph g,weight_edge t){ 
	weight_edge::iterator it;
	weight_edge::iterator it_g;
	int num;
	for(it=t.begin(); it != t.end(); it++){ //O(|V-1|)
		num=0;
	
		for(it_g=g.edges.begin(); it_g != g.edges.end(); it_g++){ //O(E)
			if(it->second.first==it_g->second.first && it->second.second == it_g->second.second){
				g.edges.erase(g.edges.begin()+num);
			}
			else{
				num++;
			}
		}
	}
	sort(g.edges.begin(),g.edges.end()); //O(ElogE)

	return g.edges;
}

vector <pair<int,vertices>> min(weight_edge g){
	int min=32767;
	int u,v;
	vector <pair<int,vertices>> :: iterator it;
	for(it = g.begin(); it!=g.end();it++){
		if(it->first<min){
			min = it->first;
			u = it->second.first;
			v = it->second.second;
		}
	}

	return {{min,{u,v}}};

}

vector <pair<int,vertices>> max(weight_edge g){
	int max = 0;
	int u,v;
	vector <pair<int,vertices>> :: iterator it;
	for(it = g.begin(); it!=g.end();it++){
		if(it->first>max){
			max = it->first;
			u = it->second.first;
			v = it->second.second;

		}
	}
	return {{max,{u,v}}};
}

vector <pair<int,vertices>> aug_BFS(vector <pair<int,vertices>> t, int u){
	int n = vertex + 1 - t.size();

	vector <pair<int,vertices>> tmp;
	int* color = new int[t.size()+n];
	int* parent = new int[t.size()+n];
	int* weight = new int[t.size()+n];

	for(int i=0 ; i < t.size()+n;i++){
		color[i]=0; // 0: white 1:gray 2:black
		parent[i]=-2;
		weight[i] = -1;
	}
	queue<int> q;
	if(color[u]==0){
		color[u]=1; //start vertex:gray;
		weight[u]=0;
		parent[u]=-1;
		q.push(u);
		while(!q.empty()){   
			int j = q.front();
			for(weight_edge::iterator it = t.begin();it!=t.end();it++){
				if(it->second.first==j && color[it->second.second]==0){
					color[it->second.second]=1;
					weight[it->second.second] = it->first;
					parent[it->second.second] = j;
					q.push(it->second.second);
					continue;
				}
				if(it->second.second==j && color[it->second.first]==0){
					color[it->second.first]==1;
					weight[it->second.first] = it->first;
					parent[it->second.first] = j;
					q.push(it->second.first);
					continue;
				}
			}
			q.pop();
			color[j] = 2;

		}
	}
	for(int index=0 ; index < t.size()+n;index++){
		if(parent[index]>0){
			if(parent[index]<index)
				tmp.insert(tmp.begin(),{weight[index],{parent[index],index}});
			else
				tmp.insert(tmp.begin(),{weight[index],{index,parent[index]}});
		}
	}	
	
	delete [] color;
	delete [] parent;
	delete [] weight;

	return tmp;

}

//O(V + E) 
vector <pair<int,vertices>> BFS(vector<pair<int,vertices>> t,int u,int v){
	int n = vertex + 1 - t.size();
	vector <pair<int,vertices>> tmp;
	int* color = new int[t.size()+n];
	int* parent = new int[t.size()+n];
	int* weight = new int[t.size()+n];

	for(int i=0 ; i < t.size()+n;i++){
		color[i]=0; // 0: white 1:gray 2:black
		parent[i]=-2;
		weight[i] = -1;
	}

	queue<int> q;
	if(color[u]==0){
		color[u]=1; //start vertex:gray;
		weight[u]=0;
		parent[u]=-1;
		q.push(u);
		while(!q.empty()){   
			int j = q.front();
			for(weight_edge::iterator it = t.begin();it!=t.end();it++){
				if(it->second.first==j && color[it->second.second]==0){
					color[it->second.second]=1;
					weight[it->second.second] = it->first;
					parent[it->second.second] = j;
					q.push(it->second.second);
				}
				if(it->second.second==j && color[it->second.first]==0){
					color[it->second.first]==1;
					weight[it->second.first] = it->first;
					parent[it->second.first] = j;
					q.push(it->second.first);
				}
			}
			q.pop();
			color[j] = 2;

		}
	}
	int index = v;
	while(parent[index]!=-1){ //worst case: O(V)
	if(parent[index]<index)
		tmp.insert(tmp.begin(),{weight[index],{parent[index],index}});
	else
		tmp.insert(tmp.begin(),{weight[index],{index,parent[index]}});
	index = parent[index];
	}
	delete [] color;
	delete [] parent;
	delete [] weight;
	return tmp;
}


vector <pair<int,vertices>> incre_min_weight_of_tree(Graph g, vector<pair<int,vertices>> t){
	int z = 32767;
	vector<pair<int,vertices>> new_t;
	vector<pair<int,vertices>>* container= new vector<pair<int,vertices>>[2] ;
	vector<pair<int,vertices>> path;
	g.edges = sort_edge(g, t);
	vector <pair<int,vertices>> :: iterator it;
	for(it = g.edges.begin(); it!=g.edges.end(); it++){  //O(|E|)
		container[0] = min(t); // |V-1|
		container[1] = max(t);
		if(it->first > container[0][0].first && it->first < container[1][0].first)
		{
			path = BFS(t,it->second.first,it->second.second); //find the path from u to v in the tree.
			t.push_back({it->first,{it->second.first,it->second.second}});
			//output3 << "(a)after adding the edge, the size of the tree: " << t.size()<< endl; //add an edge u-v. It will form a cycle C in T.
			int tmp_u;
			int tmp_v;
			int weight=32767;
			for(weight_edge :: iterator itr = path.begin(); itr != path.end(); itr++)
			{
				if(itr->first < weight)
				{
					weight = itr->first;
					tmp_u = itr->second.first;
					tmp_v = itr->second.second;
				}
			}
			if(weight >= it->first)
			{
				t.pop_back();
			}
			else{
				for(weight_edge::iterator itr = t.begin(); itr != t.end(); )
				{
					if(itr->first == weight && itr->second.first == tmp_u && itr->second.second == tmp_v){
						itr = t.erase(itr);
						break;
					}
					else{
						itr++;
					}
				}
			}
			container[0] = min(t);
			container[1] = max(t);
			if(z > container[1][0].first - container[0][0].first){
				z= container[1][0].first - container[0][0].first;
				if(z==0)
					break;
			}			
		}	
	}
	cout << "After alg2 , the size of the tree: " << t.size()<< endl;
	cout << "z value: " << z << endl;
	delete [] container;
	new_t = t;
	return new_t;
}



bool compute_all_vertex_degrees(weight_edge t,vector<int> R_prime){
	vector<list<int>> adjList(vertex+1);
	for(weight_edge::iterator it=t.begin(); it!=t.end();it++){
		adjList[it->second.first].push_back(it->second.second);
		adjList[it->second.second].push_back(it->second.first);
	}
	for(vector<int>::iterator itr = R_prime.begin(); itr != R_prime.end(); itr++){
		if(adjList[*itr].size() >=2){
			continue;
		}
		else{
			return false;
		}

	}
	return true;
}

bool adjacent(weight_edge t,int Va,int Vb){
	for(weight_edge:: iterator it = t.begin(); it != t.end(); it++){
		if(( it->second.first == Va && it->second.second == Vb) || (it->second.first == Vb && it->second.second == Va)){
			
			return false;
		}

	}
	return true;
}

weight_edge tmp;
vector<pair<int,vertices>> Hamiltonian_connected(weight_edge t,int Va , int Vb){
	int Va_prime;
	int Vb_prime;
	weight_edge Tva;
	weight_edge Tvb;
	if(t.size()!=0){
		if(adjacent(t,Va,Vb)){ 
			weight_edge P = BFS(t,Va,Vb);
			for(weight_edge::iterator it=t.begin();it!=t.end();it++){
				if(it->second.first==P[0].second.first && it->second.second == P[0].second.second){
					t.erase(it);
					break;
				}
                if(it->second.first==P[0].second.second && it->second.second == P[0].second.first){
					t.erase(it);
					break;
				}
			}
			
			Tva = aug_BFS(t,Va);
			Tvb = aug_BFS(t,Vb);
			if(Tva.size()==0){ // Tva have only one vertex and no edge.
				Va_prime = Va;
			}
			else{
				for(weight_edge::iterator it=Tva.begin(); it!=Tva.end(); it++){
					if(it->second.first==Va)
					{
						Va_prime = it->second.second;
						break;
					}
					if(it->second.second==Va)
					{
						Va_prime = it->second.first;
						break;
					}	
				}
			}
			if(P[0].second.first != Va)
				Vb_prime = P[0].second.first;
			else
				Vb_prime = P[0].second.second;

		}
		else{
			for(weight_edge::iterator it=t.begin();it!=t.end();it++){
				if(it->second.first==Va && it->second.second == Vb){
					t.erase(it);
					break;
				}
				if(it->second.first==Vb && it->second.second == Va){
					t.erase(it);
					break;
				}
			}
			Tva = aug_BFS(t,Va);
			Tvb = aug_BFS(t,Vb);
			if(Tva.size()==0){ // Tva have only one vertex and no edge.
				Va_prime = Va;
			}
			else{
				for(weight_edge::iterator it=Tva.begin(); it!=Tva.end(); it++){
					if(it->second.first==Va)
					{
						Va_prime = it->second.second;
						break;
					}
					if(it->second.second==Va)
					{
						Va_prime = it->second.first;
						break;
					}

				}
			}
			if(Tvb.size()==0){ // Tvb have only one vertex and no edge.
				Vb_prime = Vb;
			}
			else{
				for(weight_edge::iterator it=Tvb.begin(); it!=Tvb.end(); it++){
					if(it->second.first==Vb)
					{
						Vb_prime = it->second.second;
						break;
					}
					if(it->second.second==Vb)
					{
						Vb_prime = it->second.first;
						break;
					}

				}
			}
		}
		weight_edge Pva = Hamiltonian_connected(Tva,Va,Va_prime);
		weight_edge Pvb = Hamiltonian_connected(Tvb,Vb_prime,Vb);

		if(Va_prime<Vb_prime)
		{
			tmp.push_back({-1,{Va_prime,Vb_prime}});
		}
		else
		{
			tmp.push_back({-1,{Vb_prime,Va_prime}});
		}
	}
	else{
		return t;
	}
	return tmp;
}
int main(){
	int V=vertex,E=edges;
	int edges[vertex+1][vertex+1]={0};
	Graph g(V, E);
	int vertex_A=1;
	int vertex_B=2;
    ifstream file;
    file.open("case10.txt",ios::in);
	ofstream output("edge_weight.txt",ios::out);
	ofstream output1("kruskalTree.txt",ios::out);
	ofstream output2("after_alg2.txt",ios::out);
	ofstream output3("after_alg3.txt",ios::out);
    char* Covert_char;
    string  str;
    if(file.is_open()){
        cout << "open file successfully.\n";
    }
    else{

        cout<< "open file failed.\n";
    }
	vector<int> R;
	vector<int> R_prime;
	int input;
	cout << "Enter elements of the set R(1~280): (-1:exit)" << endl;
	while(input != -1){
		cin >> input;
		if(input != -1)
			R.push_back(input);
	}
	sort(R.begin(),R.end());

	cout << "Enter elements of the set R'(R' is a proper set of R): (-1:exit)" << endl; //|R\R'|>=2
	int input1;
	while(input1 != -1){
		cin >> input1;
		if(input1 != -1)
			R_prime.push_back(input1);
	}
	sort(R_prime.begin(),R_prime.end());
	vector<int> diff;
	set_difference(R.begin(),R.end(),R_prime.begin(),R_prime.end(),inserter(diff,diff.begin()));
	
    while(getline(file , str ,'\n')){

        Covert_char = strcpy(new char[str.length() + 1],str.c_str());
		int weight = atoi(Covert_char);
		if(vertex_B < vertex+1){
			output << vertex_A << "-" << vertex_B << "weight: " << weight  <<endl;
			edges[vertex_A][vertex_B] = weight;
			g.addEdge(vertex_A,vertex_B,weight);
			vertex_B++;
		}
		else{
			vertex_A++;
			vertex_B = vertex_A+1;
			output << vertex_A << "-" << vertex_B << "weight: " << weight  <<endl;
			edges[vertex_A][vertex_B] = weight;
			g.addEdge(vertex_A,vertex_B,weight);
			vertex_B++;
		}
		

    }
	vector<pair<int,vertices>> tree = g.kruskalMST();
	vector< pair<int, vertices> >::iterator it;
	for (it = tree.begin(); it != tree.end(); it++){
		output1 << it->second.first << " - " << it->second.second <<" weight:"<< it->first << endl;
	}
	
	tree = incre_min_weight_of_tree(g,tree);
	
	sort(tree.begin(),tree.end());
	for (it = tree.begin(); it != tree.end(); it++){
		output2 << it->second.first << " - " << it->second.second <<" weight:"<< it->first << endl;
	}
	
	if(compute_all_vertex_degrees(tree,R_prime)){
		cout << "Before alg3, it is a selected-internal steiner tree." << endl;
	}
	else{
		cout << "Before alg3, it is not a selected-internal steiner tree." << endl;
		random_device rd;
		default_random_engine rng(rd());
		shuffle(begin(diff),end(diff),rng);
		cout << diff[0] <<", " << diff[1] << endl;
		tree = Hamiltonian_connected(tree,diff[0],diff[1]);
		cout <<"after alg3, t.size(): " <<tree.size() << endl;
	}

	for (it = tree.begin(); it != tree.end(); it++){
		it->first = edges[it->second.first][it->second.second];
	}
	sort(tree.begin(),tree.end());
	for (it = tree.begin(); it != tree.end(); it++){
		output3 << it->second.first << " - " << it->second.second <<" weight:"<< it->first << endl;
	}
	/* Check the tree.
	if(compute_all_vertex_degrees(tree,R_prime)){
		cout << "It is a selected-internal steiner tree." << endl;
	}
	else{
		cout << "It is not a selected-internal steiner tree." << endl;
	}
	*/
	weight_edge min_value = min(tree);
	weight_edge max_value = max(tree);
	cout << "A balanced selected-internal steiner tree: " << max_value[0].first - min_value[0].first << endl;

	output3.close();
	output2.close();
	output1.close();
	output.close();
	file.close();
    return 0;
}