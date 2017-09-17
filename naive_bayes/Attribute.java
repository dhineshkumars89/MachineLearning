package assignment3;

public class Attribute {
	double[] mean;
	double[] sd;
	double[] W;
	 public Attribute(int atrnum,double min,double max,int N) 
	 {
		W = new double[N];
		mean = new double[N];
		sd = new double[N];
		double G,L=max,S=min;
		G=(L-S)/N;
		if((G<0.0001) || Double.isInfinite(G) || Double.isNaN(G))
			G=0.0001;
		//intialize all mean,sd,W
			for(int i=0;i<N;i++)
			{
				W[i] = 1.0/(double)N;
				sd[i] = 1.0;
				mean[i] = S + ((double)i * G) + (G/2);
			}
	        
	 }
	 @Override
    public String toString() {
		 StringBuilder strb = new StringBuilder();
		 for(int i=0;i<W.length;i++)
			{
				strb.append(String.format("%.4f		%.4f	%.4f",W[i],sd[i],mean[i]));
				strb.append("\n");
			}
        return strb.toString();

    }

//	Read more: http://javarevisited.blogspot.com/2012/09/override-tostring-method-java-tips-example-code.html#ixzz4Y34wNXIv
	 
	 public double[] getMean() {
			return mean;
		}
		public void setMean(double[] mean) {
			this.mean = mean;
		}
		public double[] getSd() {
			return sd;
		}
		public void setSd(double[] sd) {
			this.sd = sd;
		}
		public double[] getW() {
			return W;
		}
		public void setW(double[] w) {
			this.W = w;
		}
	 

}
