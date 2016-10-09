
public class Collision_Information {

	public double[][] Point;
	public double time;
	
	Body Object;

	Boolean Collision_Occurred;
	
	Collision_Information(int size)
	{
		Point = new double[size][3];
		time = 0;
		Collision_Occurred = false;
	}
	
	
}
