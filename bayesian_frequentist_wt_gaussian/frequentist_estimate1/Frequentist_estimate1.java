import java.util.Random;
import java.util.concurrent.TimeUnit;

public class Frequentist_estimate1 {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		long startTime = System.currentTimeMillis();
		
		Random rnd = new Random();
		 //rnd.setSeed(2000);
		 float val;
		 String str = "";
		 float numa = 0;
		 int count = 3100;
		 for(int i=0;i < count; i++)
		 {
			 val = rnd.nextFloat();
			 if(val <= 0.1)
			 {
				 str = str + 'a';
				 numa++;
			 }
			 else
			 {
				 str = str + 'b';
				 //numb++;
			 }
		 }
		 //System.out.println(str);
		 System.out.println(String.format("p(c = 'a') = %.4f",(numa/count)));
		 
		 		 
		 long endTime   = System.currentTimeMillis();
		 long totalTime = endTime - startTime;
//		 System.out.println(totalTime);
//		 System.out.println(TimeUnit.MILLISECONDS.toSeconds(totalTime));
		 //TimeUnit.MILLISECONDS.toSeconds(totalTime)
	}

}
