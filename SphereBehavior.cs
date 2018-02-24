using System.Collections;
using System.Collections.Generic;
using System.Runtime;/// <summary>
/// where is Int32 defined?
/// </summary>
using UnityEngine;

public class SphereBehavior : MonoBehaviour {

//	public ComputeShader shader;
	public float noiseSize = 0.001f;
	public float turnRate = 0.1f;
	public int paintThisEdge = 1;
	public float rateAdditionR = 0.01f;
	public float rateSubtractionG = 0.001f;
	public float rateTransform = 0.001f;
	public float rateDiffusionR = 0.001f;
	public float rateDiffusionG = 0.0001f;
	public int dropPoint = 20;
	public Color dropColor;

	private Mesh theMesh;
	private MeshFilter theMeshFilter;
	private Vector3[] vertArray;
	private float[] vertAreas;
	private Vector3[] normArray;
	private int[] trianglesArray;
	private int nV;
	private int nT;
	private int nE;

	float EdgeLength( int n1, int n2) {
		Vector3 p1 = vertArray [n1];
		Vector3 p2 = vertArray [n2];
		return Vector3.Distance (p1, p2);
	}

	float TriangleArea( int n1, int n2, int n3 ) {
		Vector3 s1 = vertArray [n2] - vertArray [n1];
		Vector3 s2 = vertArray [n3] - vertArray [n1];
		Vector3 xp = Vector3.Cross (s1, s2);
		return 0.5f * xp.magnitude;
	}

	public struct Edges {
		public int node1;
		public int node2;
		public int side1;	// index of vertex to one side of the edge, and the other
		public int side2;
		public float areaNode1;		// area of surface covered by node1 to the boundary
		public float areaNode2;
		public float boundaryLength; // so, diffusion into is proportional to surfaceLength, 
		public bool forward;
		public bool dup;
		public bool dopple;
		public Edges( int n1,int n2,int n3 ) {
			if (n1<n2) { // sort them
				node1 = n1;
				node2 = n2;
				forward = true;
			} else {
				node1 = n2;
				node2 = n1;
				forward = false;
			}
			side1 = n3;
			dup = false;
			dopple = false;
			side2 = 0;
			areaNode1 = 0;
			areaNode2 = 0;
			boundaryLength = 0;
		}

	}

	public void ComputeMyDimensions(ref Edges edg ) {
		edg.areaNode1 = TriangleArea (edg.node1, edg.side1, edg.side2);
		edg.areaNode2 = TriangleArea (edg.node2, edg.side2, edg.side1);
		edg.boundaryLength = EdgeLength (edg.side1, edg.side2);
	}


	Edges[] edges;
	void FindAllPairs() {
		int i1, i2, i3;
		int totalUnique = nT; // dec for each dup found
		Edges[] allPairs = new Edges[nT];
		for (int tn = 0; tn < nT; tn += 3) {
			i1 = trianglesArray [tn];
			i2 = trianglesArray [tn + 1];
			i3 = trianglesArray [tn + 2];
			allPairs [tn] = new Edges (i1, i2, i3);
			allPairs[tn+1] = new Edges (i2, i3, i1);
			allPairs[tn+2] = new Edges (i3, i1, i2);
		}
		for (int n = 0; n < nT-1; n++) {
			if (allPairs [n].dup)
				continue;
			for (int m = n+1; m < nT; m++) {
				if (allPairs [m].dup)
					continue;
				if ((allPairs[n].node1 == allPairs[m].node1) & (allPairs[n].node2 == allPairs[m].node2)) {
					allPairs [n].dup = true;
					allPairs [m].dup = true;
					allPairs [n].side2 = allPairs [m].side1;
					allPairs [m].side2 = allPairs [n].side1; // Not really needed, as will be excluded. 
					allPairs [m].dopple = true; // tag which one to exclude
					totalUnique -= 1;
					break;
				}
			}
		}
		nE = totalUnique;
		edges = new Edges[nE];
		int walk = 0;
		float totBound = 0;
		float totArea = 0;
		for (int j = 0; j < nT; j++) {
			if (!allPairs [j].dopple) {
				edges [walk] = allPairs [j];
				ComputeMyDimensions (ref edges[walk]);
				totBound += edges [walk].boundaryLength;
				totArea += edges [walk].areaNode1;
				vertAreas [edges [walk].node1] += edges [walk].areaNode1;
				walk += 1;
			}
		}
		Debug.Log ("n triangles " + nT);
		Debug.Log ("n edges " + nE);
		Debug.Log ("n vert " + nV);
		Debug.Log ("avg boundary length " + totBound /(1.0f* nT));
		Debug.Log ("avg area " + totArea / (1.0f*nT));
	}
		

	// 768 tri 7fps cpu i7
	// 12K verts 7680 tri 12000ms (0.1fps)
	// pretty sure edges, triangles, vertices and colors index do not match up 
	// as I am expecting.
	bool busy = false;
	private Vector4[] theColors;// x,y,z,w

	void DiffuseAlongEdgesPingPong(bool propa, int pntme ) {
		// takes r,g as time point 0, puts b,a for time point 1. eg, Ping
		// expecting the react call to react b,a and put it back in r,g. eg Pong
		Vector4 c1, c2;
		// theColors = theMesh.colors; // presume fetched and replaced outside somewhere
		int v1, v2;
		float diff;
		for (int en = 0; en < nE; en++) {
			v1 = edges [en].node1;
			v2 = edges [en].node2;
			c1 = theColors [v1];
			c2 = theColors [v2];

			// diffuse r between 1&2 at proportional rate k4
			// x,y  diffuses into z,w
			if (!propa) {
				diff = (c1.x - c2.x);// * edges[en].boundaryLength;
				c1.z -= diff * rateDiffusionR; // / vertAreas[edges[en].node1];
				c2.z += diff * rateDiffusionR; // / vertAreas[edges[en].node2];
			// and diffuse g
				diff = (c1.y - c2.y);// * edges[en].boundaryLength;
				c1.z -= diff * rateDiffusionG; // / vertAreas[edges[en].node1];
				c2.z += diff * rateDiffusionG; // / vertAreas[edges[en].node2];
			}

			theColors [v1] = c1;
			theColors [v2] = c2;
		}
	}


	void ReactBackPingPong( bool justPong = false ) {
		// expecting the react call to react z,w and put it back in x,y eg Pong
		Vector4 c1;
		float spawn;
		for (int v1 = 0; v1 < nV; v1++) {
			c1 = theColors [v1];
			c1.x = c1.z;
			c1.y = c1.w;
			if (!justPong) {
				c1.x += rateAdditionR * (1.0f - c1.x);  	// add componenet r at rate k1 
				c1.y -= rateSubtractionG * c1.y;  			// remove g at proportional rate k2
				spawn = c1.x * c1.y * c1.y * rateTransform;
				c1.y += spawn; // transform r to g at rate k3;
				c1.x -= spawn;
			}
			theColors [v1] = c1;
		}
	}

	void AssignToVerticis( ) {
		Vector3[] vs = theMesh.vertices;	
		Vector3[] ns = theMesh.normals;
		for (int v1 = 0; v1 < nV; v1++) {
			vs [v1] = ns [v1] * theColors [v1].x * noiseSize;
		}
		theMesh.vertices = vs;
	}

	void AssignRandomColors(bool solidBlue = false) {
		Color[] asgn = new Color[nV];
		for (int iC = 0; iC < nV; iC++) {
			if (solidBlue) {
				theColors [iC].x = 0.0f;
				theColors [iC].y = 0.0f;
				theColors [iC].z = 1.0f; // I lied. It is solid blue;
			} else {
				theColors [iC].x = Random.value;
				theColors [iC].y = Random.value;
				theColors [iC].z = theColors [iC].y; // 0.0f; // 1.0f * iC / nV;
				theColors [iC].w = theColors [iC].w; // 0.0f; // 1.0f * iC / nV;=
			}
			asgn [iC] = theColors [iC];
		}
		theMesh.colors = asgn;

	}

	void remesh( float threshold = 0.00001f) {
		int nVertOrg,nVertUniq;
		Vector3[] vertArray = theMesh.vertices;
		Vector3[] newVertArray;
		nVertOrg = vertArray.Length;
		int[] trianglesArray = theMesh.triangles;
		int[] newTrianglesArray;
		int[] vertIndex = new int[nVertOrg];
		int[] collapseInverseVertIndex;//



		nVertUniq = nVertOrg;

		for (int i = 0; i < nVertOrg; i++)
			vertIndex [i] = i;
		
		for (int i = 0; i < nVertOrg-1; i++) {
			if (vertIndex [i] != i)
				continue;
			for (int j = i + 1; j < nVertOrg; j++) {
				if ((vertIndex[j] == j) & (Vector3.Distance (vertArray [i], vertArray [j]) < threshold)) {
					if (Vector3.Distance (vertArray [i], vertArray [j]) != 0) { 
						Debug.Log ("D = " + Vector3.Distance (vertArray [i], vertArray [j]));
					}
					vertIndex [j] = i;
					nVertUniq -= 1;
				}
			}
		}
		// vertIndex now points into a subset of the verticis.
		// but still with indicis over the original range.
		Debug.Log ("nvertorg " + nVertOrg);
		Debug.Log ("n uniq" + nVertUniq);

		newVertArray = new Vector3[nVertUniq];
		collapseInverseVertIndex = new int[nVertUniq];
		vertAreas = new float[nVertUniq];

		int ix = 0;
		for (int i = 0; i < nVertUniq; i++) {
			vertAreas [i] = 0.0f;
			while (vertIndex [ix] != ix) { // walk up past any non-unique idx
				ix += 1; 
			}
			newVertArray [i] = vertArray [ix];
			collapseInverseVertIndex [i] = ix;
			ix += 1;
		}

		newTrianglesArray = new int[trianglesArray.Length];
		int stage1;
		for (int i = 0; i < trianglesArray.Length; i++) {
			stage1 =  vertIndex [trianglesArray [i]];
			for (int j = 0; j < nVertUniq; j++) { // aaaaaa my Brains is bleeding!
				if (collapseInverseVertIndex [j] == stage1) {
					newTrianglesArray [i] = j;
					break;
				}
			}
		}

		theMesh.Clear (); // ??
		theMesh.vertices = newVertArray;
		theMesh.triangles = newTrianglesArray;
		// Failed setting triangles. Some indices are referencing out of bounds vertices. IndexCount: 23040, VertexCount: 3842
		nV = nVertUniq;
		theMesh.RecalculateNormals ();
		theMesh.MarkDynamic ();
		theMesh.UploadMeshData (false);
		theColors = new Vector4[nV];

	}

	void Start () {
	//	nShader ();
		Color c1;
		theMeshFilter = gameObject.GetComponent<MeshFilter> ();
		theMesh = theMeshFilter.sharedMesh;
		Debug.Log ("submesh = " + theMesh.subMeshCount);
		for (int sm = 0; sm < theMesh.subMeshCount; sm++) {
			Debug.Log ("sm 0 topology is " + theMesh.GetTopology (sm));
		}
		remesh ();
		nV = theMesh.vertexCount;

//		normArray = theMesh.normals;
		vertArray = theMesh.vertices;
		trianglesArray = theMesh.triangles;
		nT = trianglesArray.Length; // why is this not seen as valid?
		FindAllPairs (); // better for shader if needed
		AssignRandomColors (false);

	}

	public void RunShader() {
	//	int kernelIndex = shader.FindKernel ("CSMain");
	//	RenderTexture tex = new RenderTexture (512, 512, 24);
	//	tex.enableRandomWrite = true;
	//	tex.Create ();

	//	shader.SetTexture (kernelIndex, "Result", tex);
	//	shader.Dispatch (kernelIndex, 512 / 8, 512 / 8, 1);

	}


	private int verti = 0;
	private Color pac = new Color (1, 0, 1, 0);
	private int seq = 0;
	void FixedUpdate () {
		transform.Rotate (new Vector3 (turnRate, 0, 0));




	}

	int cdn = 0;
	void Update () {
	//	theColors [0] = Color.white;
		theColors [dropPoint] = dropColor;
		//theColors [7] = Color.green;
		if (cdn < 0) {
			cdn = nE;
		}

		DiffuseAlongEdgesPingPong (false,cdn);
		ReactBackPingPong(true );
		AssignToVerticis ();
		cdn -= 1;
	}
}
