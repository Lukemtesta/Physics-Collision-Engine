
//Holds information about 

public class Edge_Normal_Information {

	public double[][] Normal_towards_P0;
	public double[][] Point_on_Normals_Plane;
	
	Edge_Normal_Information(int Number_of_Normals_towards_P0)
	{
		Normal_towards_P0 = new double[Number_of_Normals_towards_P0][3];
		Point_on_Normals_Plane = new double[Number_of_Normals_towards_P0][3];
	}
	
}
