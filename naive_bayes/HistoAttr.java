package assignment3;

public class HistoAttr {

	double[] bin;
	double min;
	double max;
	int N;
	 public HistoAttr(int atrnum,double min,double max,int N) 
	 {
		this.bin = new double[N];
		this.min = min;
		this.max = max;
		this.N = N;
		double G,L=max,S=min;
		G=(L-S)/N;
		//intialize all mean,sd,W
			for(int i=0;i<N;i++)
			{
				bin[i]=0;
			}
	        
	 }
	
	public double[] getBin() {
		return bin;
	}

	public void setBin(double[] bin) {
		this.bin = bin;
	}
	
	public void setidvBin(int i,double binval) {
		this.bin[i] = binval;
	}

	public double getMin() {
		return min;
	}

	public void setMin(double min) {
		this.min = min;
	}

	public double getMax() {
		return max;
	}

	public void setMax(double max) {
		this.max = max;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
