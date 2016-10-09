
//--------------------------------------------------Updated Maths 24.2.13------------------------------------------


	 /* Print Methods
	 * 
	 * Print Vector Y
	 * Print Matrix Y
	 * Print Cube Y
	 * 
	 * 
	 * Mathematical Methods:
	 * 
	 * Matrix Multiplication Y
	 * Dot Product Y
	 * Cross Product Y
	 * Magnitude Y
	 * Power Y
	 * Vector Add/Sub Y
	 * Rotation Matrix Y
	 * Binary Search Box t_collision Y
	 * Scale Matrix 
	 * Translate Matrix
	 * Box CoM Finder Y
	 * Cube Generator
	 */

	 /*Others:
	 * 
	 * Sliding
	 * Triangle Mesh
	 * Cylinder
	 */

public class Maths {
				
				
				public double[][] Quaternion_to_QuaternionMatrix(double[] q)
				{
					double a = q[0];
					double b = q[1];
					double c = q[2];
					double d = q[3];
					
					double[][] R = { {a*a + b*b - c*c - d*d, 2*b*c - 2*a*d, 2*b*d + 2*a*c},
									 {2*b*c + 2*a*d, a*a - b*b + c*c - d*d, 2*c*d - 2*a*b},
									 {2*b*d - 2*a*c, 2*c*d + 2*a*b, a*a - b*b - c*c + d*d},
									};
					
					return R;
				}
				
				public double[] Quaternion_to_Euler(double[] q)
				{
					// Attitude = Z
					// Heading = Y
					// Bank = x
					//
					//Returns Euler in Radians
					
					double qw = q[0];
					double qx = q[1];
					double qy = q[2];
					double qz = q[3];
					
					double heading = Math.atan2(2*qy*qw-2*qx*qz , 1 - (2*qy*qy) - (2*qz*qz));
					double attitude = Math.asin(2*qx*qy + 2*qz*qw);
					double bank = Math.atan2(2*qx*qw-2*qy*qz , 1 - (2*qx*qx) - (2*qz*qz));
							
					double[] Euler = {bank, heading, attitude};
					
					return Euler;
				}
				
				public double[] Matrix_Multiply(double[][] A, double[] B)
				{		
					//Initialise width = matrix B
					int w = 1;
					//System.out.println("w: " + w);
			
					//For dot multiplication, matrices must have equal height
					if(A.length != B.length)
					{
					//	System.out.println("Fail: Matrix Height Not Equal");
						return null;
					}
					int h = A.length;
					//System.out.println("h: " + h);
						
					//variable k (pointer to element along w currently being multiplied) must be equal to the width of matrix A
					int w_k = A[0].length;
					
					//Initialise return matrix to be same size as matrix B
					double[] C = new double[h];
					
					//Calculate matrix C elements based on the multiplication rule stated above
					for(int i = 0; i < h; i++)
					{
						for(int j = 0; j < w; j++)
						{
							for(int k = 0; k < w_k; k++)
							{
								C[i] += A[i][k]*B[k];
							//	System.out.println("A[" + i + "][" + k + "]: " + A[i][k] +  " B[" + k + "][" + j + "]: " + B[k][j] );
							//	System.out.println("C[i][j]: " + C[i][j]);
							}
							//System.out.println("variable printed");
						}
					}
					
					return C;
				}
				
				//Overload string is vector.R, else point.R
				
				public double[] EulerMatrix_Rotation(double[][] R, double[] p, String s)
				{
					//Multiply point p, by transform matrix R
					//
					// p' = R.p 
					//
					// p = [x y z 1] - Points
					// 
					
					double[] point = new double[4];
					point[0] = p[0]; point[1] = p[1]; point[2] = p[2];
					point[3] = 0;
					
					double[] ret = new double[3];
					ret = Matrix_Multiply(R,point);
					
					return ret;
				}
				
				public double[] EulerMatrix_Rotation(double[][] R, double[] p)
				{
					//Multiply point p, by transform matrix R
					//
					// p' = R.p 
					//
					// p = [x y z 1] - Points
					// 
					
					double[] point = new double[4];
					point[0] = p[0]; point[1] = p[1]; point[2] = p[2];
					point[3] = 1;
					
					double[] ret = new double[3];
					ret = Matrix_Multiply(R,point);
					
					return ret;
				}
				
				public double[][] Euler_to_EulerMatrix(double[] R)
				{
					//Takes in Euler's angle and creates rotation matrix
					//
					// Uses 3-2-1 body rule [ZYX] 
					//
					
					double sx = Math.sin(R[0]);
					double sy = Math.sin(R[1]);
					double sz = Math.sin(R[2]);
					double cx = Math.cos(R[0]);
					double cy = Math.cos(R[1]);
					double cz = Math.cos(R[2]);
					
					double[][] Rx = {{1, 0, 0, 0,},
									{0, cx, -sx, 0,},
									{0, sx, cx, 0},
									{0, 0, 0, 1}};
					
					double[][] Ry = {{cy, 0, -sy, 0,},
									{0, 1, 0, 0,},
									{sy, 0, cy, 0},
									{0, 0, 0, 1}};
			
					double[][] Rz =	{{cz, -sz, 0, 0,},
									{sz, cz, 0, 0,},
									{0, 0, 1, 0,},
									{0, 0, 0, 1}};
					
					return Matrix_Multiply(Rz,Matrix_Multiply(Ry,Rx));
				}
				
				public void Print_4D_Vector(String s, double[] a)
				{
					System.out.print(s + ": ");
					for(int i = 0; i < 4; i++)
					{
					System.out.print(a[i] + ", ");
					}
					System.out.println();
				}
				
				public double[] AxisAngle_Rotation(double[] AxisAngle, double[] p)
				{
					// rotates point around axisAngle
					//
					// uses rodriguez rotation rule
					//
					// p' = p.cos(angle) + (n x p).sin(angle) + n.(n.p)[1 - cos(angle)]
					//
					// n = AxisAngle Vector
					// angle = AxisAngle Rotation Coefficient
					// p = point
					// p' = rotated point
					//
					// Angle = Radians 
					
					// Derivation of Rodriguez
					// Angle Radians
					//
					// Lq(p) = p' = |p|.|q|.|q*|		-			R^4 to R^3 notation with real = 0
					//
					// q* = (q0 - q)
					//
					// Lq(p) = (|q0|^2 - |q|^2).p + 2(q.p).q + 2.q0.(q x p)
					//
					// q = cos(angle/2) + (u/|u|).sin(angle/2)
					//
					// Lq(n) = [ cos2(angle/2) - sin2(angle/2) ].p + 2.cos(angle/2).sin(angle/2).pn
					//		 = [cos(angle).p + sin(angle).pn]
					//
					// pn = p x u
					//
					// Full Proof - Quaternions & Rotations, Yan-Bin Jia, 2009
					
					
					//Protect Against incorrect entry
					if(AxisAngle.length != 4)
						return null;
					if(p.length != 3)
						return null;
					
					//Prepare Variables for Calculation
					double angle = AxisAngle[0];
					double[] p_rotated = new double[3];
					double[] n = new double[3];
					n[0] = AxisAngle[1];
					n[1] = AxisAngle[2];
					n[2] = AxisAngle[3];
					
					double[] p_dot_cos = Vector_Multiply(Math.cos(angle),p);
					double[] n_x_p = Cross_Product(n,p);
					n_x_p = Vector_Multiply(Math.sin(angle),n_x_p);
					double n_dot_p_dot_cos = Dot_Product_Right(n,p)*(1 - Math.cos(angle));
					double[] n_dot_np = Vector_Multiply(n_dot_p_dot_cos,n);
					
					p_rotated = Vector_Add(p_dot_cos,Vector_Add(n_x_p,n_dot_np));
				
					return p_rotated;
				}
				
				public double[] Euler_to_Quaternion(double[] a)
				{
					//Input Radians
					//
					// Euler Sequence: ZYX (3-2-1 Body Rule)
					
					// Attitude = Z
					// Heading = Y
					// Bank = x
					//
					
					double[] quaternion = new double[4];
			
					double heading = a[1];
					double bank = a[0];
					double attitude = a[2];
					
				    double c1 = Math.cos(heading / 2);
				    double c2 = Math.cos(attitude / 2);
				    double c3 = Math.cos(bank / 2);
				    double s1 = Math.sin(heading / 2);
				    double s2 = Math.sin(attitude / 2);
				    double s3 = Math.sin(bank / 2);
			
					quaternion[0] = (c1*c2*c3) - (s1*s2*s3);
					quaternion[1] = (s1*s2*c3) + (c1*c2*s3);
					quaternion[2] = (s1*c2*c3) + (c1*s2*s3);
					quaternion[3] = (c1*s2*c3) - (s1*c2*s3);
					
					return quaternion;
				
				}
				
				public double[] AxisAngle_to_Euler(double[] a)
				{
					//Angle radians
					// Euler Sequence: ZYX (3-2-1 Body Rule)
					
					double angle = a[0];
					double x = a[1];
					double y = a[2];
					double z = a[3];
					double x2 = x*x;
					double y2 = y*y;
					double z2 = z*z;
					
					double[] Euler = new double[3];
					
					Euler[0] = Math.atan2(y * Math.sin(angle)- x * z * (1 - Math.cos(angle)) , 1 - (y2 + z2 ) * (1 - Math.cos(angle)));
				    Euler[2] = Math.asin(x * y * (1 - Math.cos(angle)) + z * Math.sin(angle));
					Euler[1] = Math.atan2(x * Math.sin(angle)-y * z * (1 - Math.cos(angle)) , 1 - (x2 + z2) * (1 - Math.cos(angle)));
					
					return Euler;
				}
				
				public double Degrees_to_Radians(double a)
				{
					return (a*Math.PI/180);
				}
				
				public double Radians_to_Degrees(double a)
				{
					return (a*180/Math.PI);
				}
				
				public double[] Quaternion_to_AxisAngle(double[] a)
				{
					// Convert scalar Angle + Vector to Quarternion
					//
					// quaternion = cos(angle/2) + (u/|u|)sin(angle/2) 
					//
					// q(x,y,z) = q(x,y,z)/1-q(w)*q(w)
					
					double scalar;
					double[] unit_vector = new double[3];
					
					//Store Quaternion Complex numbers
					unit_vector[0] = a[1]; 
					unit_vector[1] = a[2]; 
					unit_vector[2] = a[3]; 
					
					//Convert quaternion scalar to AxisAngle rotation coefficient
					scalar = 2*Math.acos(a[0]);
					
					//map Quaternion vector to AxisAngle vector
					double s = 1/(1 - (a[0]*a[0]) );
					unit_vector = Vector_Multiply(Math.sqrt(s),unit_vector);
					
					//Create Axis Angle
					double[] AxisAngle = new double[4];
					AxisAngle[0] = scalar;
					AxisAngle[1] = unit_vector[0];
					AxisAngle[2] = unit_vector[1];
					AxisAngle[3] = unit_vector[2];
					
					return AxisAngle;
				}
				
				public double[] AxisAngle_to_Quaternion(double[] a)
				{
					// Convert AxisAngle rotation + Vector to Quaternion
					//
					// quaternion = cos(angle/2) + (u/|u|)sin(angle/2) 
					//
					// Angle is Radians
							
					//Protect a
					if(a.length != 4)
						return null;
					
					double[] quaternion = new double[4];
					double[] u = new double[3];
					
					//Store Quaternion vector as u and normalise
					u[0] = a[1];
					u[1] = a[2];
					u[2] = a[3];
					
					//double mag = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
					//u = Vector_Multiply(1/mag,u);
					
					u = Unit_Vector(u);
					
					//Convert u to Quaternion vector
					u = Vector_Multiply(Math.sin(a[0]/2),u);
					
					//System.out.println("a: " + a[0]);
					//System.out.println("a/2: " + a[0]/2);
					//System.out.println("cos(a/2): " + Math.cos(a[0]/2));
					
					//Convert AxisAngle rotation coefficient to Quaternion scalar
					quaternion[0] = Math.cos(a[0]/2);
					
					//Store quaternion
					quaternion[1] = u[0];
					quaternion[2] = u[1];
					quaternion[3] = u[2];
					
					return quaternion;
				}
				
				public double[] Quaternion_Multiplication(double[] p, double[] q)
				{
					// Quarternion: Scalar + Vector
					//
					// Input: Angle + vector of rotation
					//
					// p vector = p'
					// p = p0 + p' 
					//
					// pq = sequence of q then p rotation
					//
					
					double[] q_vector = new double[3];
					double[] p_vector = new double[3];
					double result_scalar;
					double[] result_vector = new double[3];
					double[] result = new double[4];
					
					//Protection against invalid quarternion entries
					if( (q.length == p.length) && (q.length == 4) )
					{
						q_vector[0] = q[1];
						q_vector[1] = q[2];
						q_vector[2] = q[3];
						p_vector[0] = p[1];
						p_vector[1] = p[2];
						p_vector[2] = p[3];
					}
					else
						return null;
					
					result_scalar = (q[0]*p[0]) - Dot_Product_Right(p_vector,q_vector);
			
					//System.out.println("p.q: " + Dot_Product_Right(p_vector,q_vector));
					//System.out.println("p0.q0: " + (q[0]*p[0]));
					//System.out.println("scalar: " + result_scalar);
					
					double[] p0_x_q = Vector_Multiply(q[0],p_vector);
					double[] q0_x_p = Vector_Multiply(p[0],q_vector);
					double[] p_x_q = Cross_Product(p_vector, q_vector);
					
					result_vector = Vector_Add(Vector_Add(p0_x_q, q0_x_p),p_x_q);
					
					result[0] = result_scalar;
					result[1] = result_vector[0];
					result[2] = result_vector[1];
					result[3] = result_vector[2];
					
					return result;
					
				}
				
				//--------------------------------------END OF NEW METHODS---------------------------------------------
				
				//Returns normal to face
				//
				//Inputs: Face involved in calculation [vertex][co-ordinates]
				//
				//Returns Normal to plane.
				
				public static double[] Normal_to_Plane(double[][] Box)
				{
					int j = 1; 
					double[] a = new double[3]; 
					double[] b = new double[3];
				//	a[0] = b[0] = a[1] = b[1] = a[2] = b[2] = 1;
					
					if(Box.length < 2)
						return null;
			
						a = Vector_Subtract(Box[0],Box[j]);				// Take edge 0 as reference vertex.
						b = Vector_Subtract(Box[0],Box[j+1]);
			
					return Cross_Product(a,b);
				}
			
				//Sign magnitude returns the sign (as normalised scalar) of an input
				
				public static double Sign_Magnitude_Scalar(double a)
				{
					a = NaN_Checker(a / Magnitude_Scalar(a));
					return a;
				}
				
				//Takes input radius = length of cube edge
				// Overload = Translation of centre included in cube generation
				//
				// Outputs normalised cube at origin as 24 vertices
				// Cube[face No][Vertex No][plane]
				
				static public double[][][] Cube_Generation_3D(double radius)
				{
					double[] T = {0,0,0};
					return Cube_Generation_3D(radius, T);
				}
			
				static public double[][][] Cube_Generation_3D(double radius, double[] T)
					{
						//Define Face Matrix Defining the combination of vertex number against the face the vertex appears in
						int[][] Face_Matrix = 	{		    {0, 1, 5, 4},
															{3, 2, 6, 7},
															{1, 3, 7, 5},
															{2, 0, 4, 6},
															{0, 2, 3, 1},
															{4, 5, 7, 6},
												};
						
						
								//Define Centre of Mass 
								double[] CoM = new double[3];
								CoM[0] = 0; CoM[1] = 0; CoM[2] = 0;
								
								//Define Binary List of all Cube vertices x8	[vertex][plane]
								double[][] Vn = new double[8][3];
								
								//Define Cube [face][vertex number][Plane]
								double[][][] Cube = new double[6][4][3];
								
								//Generate binary style combination of cube vertices, Vn at origin
								//
								// & = bitwise AND operator
								//
								// Vn = [ r.(-1)^n&1, r.(-1)^n&2, r.(-1)^n&4 ]
								//
								for(int n = 0; n <= 7; n++)
								{
									Vn[n][0] = T[0] + radius*Power((-1),(n&1));
									Vn[n][1] = T[1] + radius*Power((-1),(n&2));
									Vn[n][2] = T[2] + radius*Power((-1),(n&4));
								}
								
								// Use Face_Matrix to determine vertices defining each face
								// face x6 vertices x4
								for(int i = 0; i < 6; i ++)
								{
									for(int j = 0; j < 4; j++)
									{
										Cube[i][j][0] = Vn[Face_Matrix[i][j]][0];
										Cube[i][j][1] = Vn[Face_Matrix[i][j]][1];
										Cube[i][j][2] = Vn[Face_Matrix[i][j]][2];
									}
								}
								
								return Cube;
					}
			
				public double[] CoM_Box_3D(double[][][] Box)
					{
									// Calculate Cube's CoM
									//
									// CoM = [ (Ax - Bx), (Ay - Cy), (Az - Dz) ] 
									//
									// A,B,C,D represent vertices of the box where same plane components are not equal (eg. x of A != x of B)
									//
									// Maximum of 2 faces are required to get change in x, y and z. This holds true for all signs of the two vertex signs.
									//
									double[] CoM = new double[3];
									for(int j = 0; j < 3; j++)
									{
										if(Box[0][0][0] != Box[1][j][0])
											CoM[0] = Box[0][0][0] - ((Box[0][0][0] - Box[1][j][0]) / 2);
										if(Box[0][0][1] != Box[1][j][1])
											CoM[1] = Box[0][0][1] - ((Box[0][0][1] - Box[1][j][1]) / 2);
										if(Box[0][0][2] != Box[1][j][2])
											CoM[2] = Box[0][0][2] - ((Box[0][0][2] - Box[1][j][2]) / 2);
									}
									
									return CoM;
					}
			
				// Dot multiplies any two size matrices. If one matrix, A, is larger than other, the submatrix of A (size of B) is calculated
				//
				// C[i][j] += A[i][k]*B[k][j]
				//
				// Computation: O(n^3)
				//
				// Assumes A and B are of equal height.
				// Result matrix M is size of matrix B
				//
				// Return: Returns C = null if matrix are not equal n x n sizes
				
				public double[][] Matrix_Multiply(double[][] A, double[][] B)
				{		
					//Initialise width = matrix B
					int w = B[0].length;
					//System.out.println("w: " + w);
			
					//For dot multiplication, matrices must have equal height
					if(A.length != B.length)
					{
					//	System.out.println("Fail: Matrix Height Not Equal");
						return null;
					}
					int h = A.length;
					//System.out.println("h: " + h);
						
					//variable k (pointer to element along w currently being multiplied) must be equal to the width of matrix A
					int w_k = A[0].length;
					
					//Initialise return matrix to be same size as matrix B
					double[][] C = new double[h][w];
					
					//Calculate matrix C elements based on the multiplication rule stated above
					for(int i = 0; i < h; i++)
					{
						for(int j = 0; j < w; j++)
						{
							for(int k = 0; k < w_k; k++)
							{
								C[i][j] += A[i][k]*B[k][j];
							//	System.out.println("A[" + i + "][" + k + "]: " + A[i][k] +  " B[" + k + "][" + j + "]: " + B[k][j] );
							//	System.out.println("C[i][j]: " + C[i][j]);
							}
							//System.out.println("variable printed");
						}
					}
					
					return C;
				}
			
				//Prints Matrix to screen
				// Overload: Name of matrix
				
				public void Print_Matrix(double[][] M)
				{
					for(int i = 0; i < M.length; i ++)
					{
						System.out.print("[");
						for(int j = 0; j < M[0].length; j++)
						{
							System.out.print(M[i][j] + ", ");
						}
						System.out.println("]");
					}		
				}
			
				public void Print_Matrix(String s, double[][] M)
				{
					System.out.println(s);
					for(int i = 0; i < M.length; i ++)
					{
						System.out.print("[");
						for(int j = 0; j < M[0].length; j++)
						{
							System.out.print(M[i][j] + ", ");
						}
						System.out.println("]");
					}		
				}
				
				// Scales cube's current vertices by S
				//
				// Outputs translated cube.
				// Cube[face No][Vertex No][plane]
				
				public double[][][] Scale_Box_3D(double[][][] Cube, double[] S)
				{
					for(int i = 0; i < 6; i ++)
					{
						for(int j = 0; j < 4; j++)
						{
							Cube[i][j][0] = Cube[i][j][0] * S[0];
							Cube[i][j][1] = Cube[i][j][1] * S[1];
							Cube[i][j][2] = Cube[i][j][2] * S[2];
						}
					}
					
					return Cube;
				}
			
					
				public static double Power(double a, int b)
				{
					if(b == 0)
						a = (a/Mag(a));
					else if(b == 1)
						a = 1;
					else if(b > 1)
					{
						double ret = 1;
						for(int i = 1; i<= b; i++)
							ret = ret * a;
						a = ret;
					}
						
					return a;
				}
				
				public static double Mag(double a)
				{
					return a*a;
				}
								
				// Print's cube's vertices to screen. Face per line.
				//
				// Overload: Input String is name of cube being printed (optional).
				//
				// Cube[face No][Vertex No][plane]
				
				public void Print_2D_Box_Details(double[][] C)
				{
					Print_2D_Box_Details(C, "");
				}
				
				public void Print_2D_Box_Details(double[][] C, String s)
				{
					System.out.println(s);
			
						for(int j = 0; j < 4; j++)
						{
							System.out.print("(" + C[j][0] + "," + C[j][1] + "," + C[j][2] + ") ");
						}
						System.out.println();
				}
				
				public void Print_3D_Box_Details(double[][][] C)
				{
					Print_3D_Box_Details(C,"");
				}
				
				public void Print_3D_Box_Details(double[][][] C, String s)
				{
					System.out.println(s);
					
					for(int i = 0; i < 6; i++)
					{
						for(int j = 0; j < 4; j++)
						{
							System.out.print("(" + C[i][j][0] + "," + C[i][j][1] + "," + C[i][j][2] + ") ");
						}
						System.out.println();
					}
				}
				
				//Returns whether the two vectors are identical. Vectors must be same size
				
				public static Boolean Vector_Sign_Compare(double[] A, double[] B)
				{
					int count = 0;
					
					Print_Vector("A",A);
					Print_Vector("B",B);
					
					for(int i = 0; i < A.length; i++)
					{
						if(Sign_Magnitude_Scalar(A[i]) == Sign_Magnitude_Scalar(B[i]))
							count++;
					}
					
					if (count == A.length)
					{
						System.out.println("returning true");
						return true;
					}
					
					System.out.println("A.length: " + A.length);
					System.out.println("count: " + count);
						
					if(count == A.length)
						return true;
					else
						return false;
				}
		
		
				//Returns the distance between plane and parallel infinite plane passing through origin
				//
				// Plane Normal*[x y z] = D
				//
				// return |D|.N
				//
				
				public double[] D_Plane(double[] Normal, double[] Point)
				{
					double D = (Normal[0]*Point[0]) + (Normal[1]*Point[1]) + (Normal[2]*Point[2]);
					
					return Vector_Multiply(D,Normal);
				}
				
				//Returns whether Point passes through plane - True or False
				//
				// Find whether point passes through plane edge
				//
				// P(0) = P(0) - Nedge.Dedge;
				// P(1) = P(1) - Nedge.Dedge
				//
				// [P(0)/|P(0)|].Nedge != [P(1)/|P(1)|].Nedge			-				Point Passes through Plane!!!!
				//
				// If passes through a plane edge, then collision may be possible.
				//
				
				static Boolean Point_Pass_Through_Plane(double[] Plane_Normal, double[] Distance_Plane_to_Origin, double[] P0, double[] P1)
				{
					Boolean Pass_Through_Infinite_Plane = false;
					
					double[] P0_act = Vector_Subtract(P0,Distance_Plane_to_Origin);
					double[] P1_act = Vector_Subtract(P1,Distance_Plane_to_Origin);

					//Check if sign t = 0 is same as sign t = 1
					if(Vector_Sign_Compare(P0_act, P1_act) == false)
					{
						Print_Vector("Passes Through N",Plane_Normal);
						Pass_Through_Infinite_Plane = true;
					}

					return Pass_Through_Infinite_Plane;
				}
				
				
				static double[] D_Plane_to_Parallel_Plane_Through_Origin(double[] Normal, double[] Point)
				{
					double D = Magnitude_Scalar((Normal[0]*Point[0]) + (Normal[1]*Point[1]) + (Normal[2]*Point[2]));
					
					return Vector_Multiply(D,Normal);
				}
				
		
		
				//Returns: Array Corresponding to whether each element is pointing towards P0
				//
				//
				//Find which edge normal is pointing towards P0 - must be edge normal P0 intersects with
				//May be more than 1 normal. In this case, repeat below steps for each edge normal (2 max for 2D planes)
				//
				// (P0 - Nedge.Dedge).Nedge >= 0			-			Pointing towards point
				//
				// Dedge = Nedge*[x,y,z]
				//
				
	
				public static Boolean[] Edge_Normals_Pointing_Towards_P0(double[][] Plane, double[][] Plane_Edge_Normals, double[] P0)
				{
					//Find which edge normal is pointing towards P0 - must be edge normal P0 intersects with
					//May be more than 1 normal. In this case, repeat below steps for each edge normal (2 max for 2D planes)
					//
					// (P0 - Nedge.Dedge).Nedge >= 0			-			Pointing towards point
					//
					// Dedge = Nedge*[x,y,z]
					//
					double[] Dedge = new double[Plane.length];	
					
					//Initialise normal count obeying condition as -1. 
					Boolean[] normal_towards_P0 = new Boolean[Plane.length];
					int counter = 0;
					
					for(int i = 0; i < Plane.length; i++)
					{
						normal_towards_P0[i] = false;
					}
					
					int j = 0;
					for (int i = 0; i < Plane.length; i++)
					{
						j = i + 1;
						if(j == Plane.length)
							j = 0;
						
						Dedge[i] = (Plane_Edge_Normals[i][0]*Plane[j][0]) + (Plane_Edge_Normals[i][1]*Plane[j][1]) + (Plane_Edge_Normals[i][2]*Plane[j][2]);
						double[] Nedge_dot_Dedge = Vector_Multiply(Dedge[i],Plane_Edge_Normals[i]);
						double[] P0_minus_Dedge = Vector_Subtract(P0,Nedge_dot_Dedge);
						
						//Print_Vector("Nedge.Dedge",Nedge_dot_Dedge);
						
						//Perform Check
						//
						//(P0 - Nedge.Dedge).Nedge >= 0
						double[] test = Dot_Product(P0_minus_Dedge,Plane_Edge_Normals[i]);
						
					//	Print_Vector("P0 - N.Dedge",P0_minus_Dedge);
					//	Print_Vector("test",test);
						
						//As Edge Normals are pointing towards CoM of plane, Edge normal is pointing in opposite direction to P0. So...
						//
						// Condition <= 0
						//
						// We would like normal of edge to be pointing towards P0. Hwoever the edge normal is pointing towards
						// The plane CoM. If assuming no collision at t = 0, then the normal of the edge closest to P0 will always
						// be initially pointing away from P0. 
						//
						//
						// As could be more than 1 normal, count number of normals and their position in array for later storage
						if( test[0] <= 0 && test[1] <= 0 && test[2] <= 0 )
						{
							counter++;
							normal_towards_P0[i] = true;
						}
						else
							normal_towards_P0[i] = false;
					}
					
					return normal_towards_P0;
				}
		
					//Calculates Normal to faces of a plane
					//
					//Return: double[edge No][co-ordinate] Array of normals orthogonal to plane edges. No of normals = No plane edges
					//
					//Inputs: double[edge No][co-ordinate] Plane, double[] Normal to Plane
					
					public static double[][] Calculate_Normal_to_Plane_Edges(double[][] Plane, double[] Normal_to_Plane)
					{
						double[][] Finite_Plane_Normals = new double[Plane.length][3];
						for(int i = 0; i <  Plane.length; i++)
						{
							int j = i + 1;
							if(j == Plane.length)
							{
						//		System.out.println("i: " + i);
								j = 0;
							}
							//System.out.println(Finite_Plane_Normals.length);
							Finite_Plane_Normals[i] = Unit_Vector(Cross_Product(Normal_to_Plane, Vector_Subtract(Plane[j],Plane[i])));
							Print_Vector(("N" + i),Finite_Plane_Normals[i]);
						}
						
						return Finite_Plane_Normals;
					}
					
					//Calculates Cross Product of two 3D Vectors
					//
					//Input: A 3D vector (x,y,z) as double
					//		 B 3D vector (x,y,z) as double
					//
					//Return: C a 3D Vector (AxB)
					//
					// AxB = (AyBz - AzBy)i - (AxBz - AzBx)j + (AxBy - AyBx)k
					
				public static double[] Cross_Product(double[] A, double[] B)
					{
						double[] AxB = new double[3];
						
						AxB[0] = (A[1]*B[2]) - (A[2]*B[1]);
						AxB[1] = -((A[0]*B[2]) - (A[2]*B[0]));
						AxB[2] = (A[0]*B[1]) - (A[1]*B[0]);
						
						return AxB;
					}
					
				public static double NaN_Checker(double N)
					{
						if( Double.isNaN(N) == true)
							N = 0;
						
						return N;
					}
				
					
				//Calculates Two 3D Vectors Added
				//
				//Input: A 3D vector (x,y,z) as double
				//		 B 3D vector (x,y,z) as double
				//
				//Return: C a 3D Vector (A-B)
				
				public static double[] Vector_Add(double[] S, double[] A)
				{
					double[] C = new double[3];
					
					C[0] = A[0] + S[0];
					C[1] = A[1] + S[1];
					C[2] = A[2] + S[2];
					
					return C;
				}
				
				public double[] Vector_Add_Scalar(double[] S, double A)
				{
					double[] C = new double[3];
					
					C[0] = A + S[0];
					C[1] = A + S[1];
					C[2] = A + S[2];
					
					return C;
				}
				
				public double[] Dot_Product_Divide(double[] S, double[] A)
				{
					double[] C = new double[3];
					
					C[0] = A[0] / S[0];
					C[1] = A[1] / S[1];
					C[2] = A[2] / S[2];
					
					return C;
				}
				
				//Calculates 3D Vector Multiplied by a scalar
					//
					//Input: A scalar
					//		 B 3D vector (x,y,z) as double
					//
					//Return: C a 3D Vector (A-B)
					
				public static double[] Vector_Multiply(double S, double[] A)
					{
						double[] C = new double[3];
						
						C[0] = A[0]*S;
						C[1] = A[1]*S;
						C[2] = A[2]*S;
						
						return C;
					}
				
				public double[] Vector_Divide(double S, double[] A)
				{
					double[] C = new double[3];
					
					C[0] = A[0]/S;
					C[1] = A[1]/S;
					C[2] = A[2]/S;
					
					return C;
				}
					
				
				//Method prints vector to screen
				//
				//Inputs:	 Name of Vector
				//			 Vector as 3D double (x,y,z)
				//
				//No Return
				
				public static void Print_Vector(String name, double[] L)
				{
					System.out.print(name + ": ");
					System.out.print(L[0] + ", ");
					System.out.print(L[1] + ", ");
					System.out.println(L[2]);
				}
				
				//Calculates Unit Vector of 3D Vector
					//
					//Input:  A 3D vector (x,y,z) as double
					//
					//Return: B a 3D Vector (Unit Vector of A)
					//
					// 		  B = A/|A|
					
				public static double[] Unit_Vector(double[] A)
					{
						double[] C = new double[3];
						double Mag = Magnitude_3D(A);
						
						// Includes NaN Protection
						
						if(A[0] != 0)
							C[0] = A[0]/Mag;
						else
							C[0] = 0;
						
						if(A[1] != 0)
							C[1] = A[1]/Mag;
						else
							C[1] = 0;
						
						if(A[2] != 0)
							C[2] = A[2]/Mag;
						else
							C[2] = 0;
						
						return C;
					}
				
				//Calculates Two 3D Vectors Subtraction
				//
				//Input: A 3D vector (x,y,z) as double
				//		 B 3D vector (x,y,z) as double
				//
				//Return: C a 3D Vector (A-B)
				
				public static double[] Vector_Subtract(double[] A, double[] B)
				{
					double[] C = new double[A.length];
					
					for(int i = 0; i < C.length; i++)
					{
						C[i] = A[i] - B[i];
					}
					
					return C;
				}
				
				//Calculates Magnitude_3D of a 3D Vector
				//
				//Inputs:	A = 3D vector as double x,y,z
				//
				//Returns:  Magnitude_3D as double
				//
				
				public static double Magnitude_3D (double[] A)
				{
					// Mag = |A| = sqrt(x*x + y*y + z*z)
					
					double Mag = Math.sqrt((A[0]*A[0]) + (A[1]*A[1]) + (A[2]*A[2]));
					
					return Mag;
				}
				
				public static double Magnitude_Scalar (double A)
				{
					double Mag = Math.sqrt(A*A);
					
					return Mag;
				}
				
				//Calculates Power of scalar
				//Indici must be a positive integer (Cannot Compute Decimal Indicies)
				//
				//Inputs:	A = Scalar of any size
				//			B = multiplier
				
				public double Power_Scalar (double A, int B)
				{
					// Power_Scalar  = 0  => 1
					// Power_Scalar > 0 => A *= A*A looped Power_Scalar times
					// Power_Scalar < 0 => 1 / (A *= A*A looped Power_Scalar times)
					
					double Multiplier = A;
					
					if(B>0)			// Positive Power_Scalar => Normal Power_Scalar
					{
						for(int i = 1; i < B; i++)
							A *= Multiplier;
					}
					
					else if(B<0)	// Negative Power_Scalar => 1/Normal Power_Scalar
					{
						B = (int) Math.sqrt(B*B);
						
						for(int i = 1; i < B; i++)
							A *= Multiplier;
						A = (1 / A);
					}
					else			// Power_Scalar of 0 = 1
					{	
						if(A < 0)
							A = -1;
						else if (A > 0)
							A = 1;
						else
							A = 0;
					}
					
					return A;
				}
				
				
				//Calculates dot product
				//
				//Inputs:	A = Vector 1 as double x,y,z			-> 3 element array
				//			B = Vector 2 as double x,y,z			-> 3 element array
				//
				//Returns:  C = Vectors [x1*x2, y1*y2, z1*z2]
				//
				public static double[] Dot_Product (double[] A, double[] B)
				{
					//vectors [a,b ,c] and [d, e, f]
					//C = A.B
					//C = [a*d, b*e, c*f]
					
					double[] C = new double[3];
					C[0] = A[0]*B[0];
					C[1] = A[1]*B[1];
					C[2] = A[2]*B[2];
					
					return C;
				}
				
				public double Dot_Product_Right (double[] A, double[] B)
				{
					//vectors [a,b ,c] and [d, e, f]
					//C = A.B
					//C = [a*d, b*e, c*f]
					
					double[] C = new double[3];
					C[0] = A[0]*B[0];
					C[1] = A[1]*B[1];
					C[2] = A[2]*B[2];
					
					double val = C[0] + C[1] + C[2];
					
					return val;
				}
		
}
