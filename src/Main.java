import com.vm.jcomplex.Complex;
import javax.sound.sampled.*;
import javax.swing.*;
import java.awt.event.*;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class Main
{
	public static void main(String[] args)
	{
		Analyzer analyzer = new Analyzer();
		analyzer.LoadGUI();
	}
}

class Analyzer
{
	Boolean running = true;
	Boolean searching = true;
	Boolean fast = false;
	Boolean hamming = false;
	Boolean sub = false;
	MyFrame myFrame = new MyFrame();

	private AudioFormat getFormat()
	{
		float sampleRate = 44100;
		int sampleSizeInBits = 8;
		int channels = 1;          //Монофонический звук
		boolean signed = true;     //Флаг указывает на то, используются ли числа со знаком или без
		boolean bigEndian = true;  //Флаг указывает на то, следует ли использовать обратный (big-endian) или прямой (little-endian) порядок байтов
		return new AudioFormat(sampleRate, sampleSizeInBits, channels, signed, bigEndian);
	}

	public void LoadGUI()
	{
		JButton searchBut = new JButton("Search");
		searchBut.setBounds(202,100,95,30);
		myFrame.add(searchBut);
		JButton recalculateBut = new JButton("Recalculate");
		recalculateBut.setBounds(202,200,95,30);
		myFrame.add(recalculateBut);
		searchBut.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e)
			{
				searching = true;
				myFrame.removeAll();
				JButton hammingBut = new JButton("Hamming search");
				hammingBut.setBounds(100,100,150,30);
				myFrame.add(hammingBut);
				JButton subBut = new JButton("Subtraction search");
				subBut.setBounds(100,150,150,30);
				myFrame.add(subBut);
				JButton fastBut = new JButton("Fast search");
				fastBut.setBounds(100,200,150,30);
				myFrame.add(fastBut);
				searching = true;
				fast = false;
				hamming = false;
				hammingBut.addActionListener(new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						hamming = true;
						myFrame.removeAll();
						Analyze();
					}
				});
				subBut.addActionListener(new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						sub = true;
						myFrame.removeAll();
						Analyze();
					}
				});
				fastBut.addActionListener(new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						fast = true;
						myFrame.removeAll();
						Analyze();
					}
				});
			}
		});
		recalculateBut.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e)
			{
				searching = false;
				Analyze();
			}
		});
	}

	public void Analyze()
	{
		try
		{
			ArrayList<String> hashes = new ArrayList<>();
			ArrayList<String> freqs = new ArrayList<>();

			if (searching)
			{
				JButton b=new JButton("Stop");
				b.setBounds(202,100,95,30);
				myFrame.add(b);
				b.addActionListener(new ActionListener()
				{
					public void actionPerformed(ActionEvent e)
					{
						running = false;
					}
				});

				byte[] buffer = new byte[2048];
				final AudioFormat format = getFormat(); //Заполнить объект класса AudioFormat параметрами
				DataLine.Info info = new DataLine.Info(TargetDataLine.class, format);
				final TargetDataLine line = (TargetDataLine) AudioSystem.getLine(info);
				line.open(format, buffer.length);
				line.start();

				ByteArrayOutputStream out = new ByteArrayOutputStream();
				running = true;
				try
				{
					while (running)
					{
						int count = line.read(buffer, 0, buffer.length);
						if (count > 0)
						{
							out.write(buffer, 0, count);
						}
						Complex[][] results = Transform(out);

						Determinator determinator = new Determinator();
						ArrayList<String>[] determinatedData = determinator.Determinate(results);
						hashes = determinatedData[0];
						freqs = determinatedData[1];
					}
					out.close();
					search(hashes, freqs);
				}
				catch (IOException e)
				{
					System.err.println("I/O problems: " + e);
					System.exit(-1);
				}
			}
			else
			{
				final AudioFormat format = getFormat();
				File musicBase = new File(".\\Music");
				String[] music = musicBase.list();
				for (String m : music)
				{
					Path path = Paths.get(".\\Music\\" + m);
					File file = new File(String.valueOf(path));
					AudioInputStream in = AudioSystem.getAudioInputStream(file);
					AudioInputStream convert = AudioSystem.getAudioInputStream(format, in);
					byte[] data = new byte[convert.available()];

					ByteArrayOutputStream out = new ByteArrayOutputStream();
					try
					{
						int count = convert.read(data, 0, data.length);
						if (count > 0)
						{
							out.write(data, 0, count);
						}
						Complex[][] results = Transform(out);

						Determinator determinator = new Determinator();
						ArrayList<String>[] determinatedData = determinator.Determinate(results);
						hashes = determinatedData[0];
						freqs = determinatedData[1];
						out.close();
						while (hashes.get(hashes.size() - 1).equals("00000000000"))
						{
							hashes.remove(hashes.size() - 1);
						}
						while (hashes.get(0).equals("00000000000"))
						{
							hashes.remove(0);
						}
						path = Paths.get(".\\HashDB\\" + m).normalize();
						String st = path.toString();
						int index = st.indexOf(".");
						String name = st.substring(0, index);
						System.out.println(name);
						path = Paths.get(name + ".txt");
						Files.write(path, hashes, StandardCharsets.UTF_8);

						while (freqs.get(freqs.size() - 1).equals("0 0 0 0 0"))
						{
							freqs.remove(freqs.size() - 1);
						}
						while (freqs.get(0).equals("0 0 0 0 0"))
						{
							freqs.remove(0);
						}
						path = Paths.get(".\\DB\\" + m).normalize();
						index = path.toString().indexOf(".");
						name = path.toString().substring(0, index);
						System.out.println(name);
						path = Paths.get(name + ".txt");
						Files.write(path, freqs, StandardCharsets.UTF_8);
					}
					catch (IOException e)
					{
						System.err.println("I/O problems: " + e);
						System.exit(-1);
					}
				}

			}
		}
		catch (Exception e)
		{
			System.err.println(e + " in Analyze");
			System.exit(-1);
		}
	}

	private void search(ArrayList<String> hashes, ArrayList<String> freqs) throws IOException
	{
		File file = new File(".\\DB");
		String[] db = file.list();

		//Поиск совпадения
		int foundHHash = 0;
		int distanceHHash = 10000000;

		int foundSFr = 0;
		long distanceSFr = Long.MAX_VALUE;

		int foundOf = 0;
		int matches = 0;

		for(int i = 0; i < Objects.requireNonNull(db).length; ++i)
		{
			////////////////////////////////////////////////////////// Поиск по хешам
			Path pathHash = Paths.get(".\\HashDB\\" + db[i]);
			List<String> readHash = Files.readAllLines(pathHash);

			Path pathFr = Paths.get(".\\DB\\" + db[i]);
			List<String> readFr = Files.readAllLines(pathFr);
			if (hamming) {
				if (Files.lines(pathHash).count() > hashes.size()) {
					for (int j = 0; j < Files.lines(pathHash).count() - hashes.size(); ++j) {
						int d = 0;
						for (int k = 0; k < hashes.size(); ++k) {
							d += HammingDistance.countDistance(Long.parseLong(readHash.get(k + j)), Long.parseLong(hashes.get(k)));
						}
						if (d < distanceHHash) {
							foundHHash = i;
							distanceHHash = d;
						}
						d = 0;
					}
				} else {
					for (int j = 0; j < hashes.size() - Files.lines(pathHash).count(); ++j) {
						int dH = 0;
						for (int k = 0; k < readHash.size(); ++k) {
							dH += HammingDistance.countDistance(Long.parseLong(readHash.get(k)), Long.parseLong(hashes.get(k + j)));
						}
						if (dH < distanceHHash) {
							foundHHash = i;
							distanceHHash = dH;
						}
					}
				}
			}
			////////////////////////////////////////////////////////// Поиск по частотам
			if (sub) {
				if (Files.lines(pathFr).count() > freqs.size()) {
					for (int j = 0; j < Files.lines(pathFr).count() - freqs.size(); ++j) {
						int dS = 0;
						for (int f = 0; f < freqs.size(); ++f) {
							String[] wordsF = freqs.get(f).split("\\s+");
							String[] wordsR = readFr.get(f + j).split("\\s+");
							for (int k = 0; k < 5; ++k) {
								dS += SubtractionDistance.countDistance(Long.parseLong(wordsF[k]), Long.parseLong(wordsR[k]));
							}
						}
						if (dS < distanceSFr) {
							foundSFr = i;
							distanceSFr = dS;
						}
					}
				} else {
					for (int j = 0; j < freqs.size() - Files.lines(pathFr).count(); ++j) {
						int dS = 0;
						for (int f = 0; f < readFr.size(); ++f) {
							String[] wordsF = freqs.get(f + j).split("\\s+");
							String[] wordsR = readFr.get(f).split("\\s+");
							for (int k = 0; k < 5; ++k) {
								dS += SubtractionDistance.countDistance(Long.parseLong(wordsF[k]), Long.parseLong(wordsR[k]));
							}
						}
						if (dS < distanceSFr) {
							foundSFr = i;
							distanceSFr = dS;
						}
					}
				}
			}
			if (fast) {
				////////////////////////////////////////////////////////// Быстрый поиск по смещению
				ArrayList<Integer> offset = new ArrayList<>();
				for (int j = 0; j < readHash.size(); ++j) {
					for (int k = 0; k < hashes.size() && k <= j; ++k) {
						if (readHash.get(j).equals(hashes.get(k))) {
							offset.add(j - k);
						}
					}
				}
				for (int j = 0; j < offset.size() - 1; ++j) {
					int match = 0;
					for (int k = j + 1; k < offset.size(); k++) {
						if (offset.get(j).equals(offset.get(k))) {
							++match;
						}
					}
					if (match > matches) {
						matches = match;
						foundOf = i;
					}
				}
			}
		}
		System.out.println("Hash Hamming: " + distanceHHash + " " + foundHHash);
		System.out.println("Frequencies Sub: " + distanceSFr + " " + foundSFr);
		System.out.println("Offset: " + matches + " " + db[foundOf]);
	}

	private Complex[][] Transform(ByteArrayOutputStream out)
	{
		byte[] audio = out.toByteArray();
		final int totalSize = audio.length;
		int amountPossible = totalSize/ AnalyzeData.CHUNK_SIZE;
		Complex[][] results = new Complex[amountPossible][];

		//Для всех кусков:
		for(int i = 0; i < amountPossible; i++)
		{
			Complex[] complex = new Complex[AnalyzeData.CHUNK_SIZE];
			for(int j = 0; j < AnalyzeData.CHUNK_SIZE; j++)
			{
				complex[j] = new Complex(audio[(i * AnalyzeData.CHUNK_SIZE) + j], 0);
			}
			//Быстрое преобразование фурье
			results[i] = FFT.fft(complex);
		}
		return results;
	}
}

class FFT
{
	public static Complex[] fft(Complex[] x)
	{
		int n = x.length;
		// Базовый случай
		if (n == 1) return new Complex[] { x[0] };

		// проверка n - степень 2, для алгоритма Кули — Тьюки
		if (n % 2 != 0)
		{
			throw new IllegalArgumentException("n не делится на 2");
		}

		// FFT для четных
		Complex[] half = new Complex[n/2];
		for (int k = 0; k < n/2; k++)
		{
			half[k] = x[2*k];
		}
		Complex[] evenFFT = fft(half);

		// повтор FFT для нечетных
		// Повторо работаем с массивом half
		for (int k = 0; k < n/2; k++)
		{
			half[k] = x[2*k + 1];
		}
		Complex[] oddFFT = fft(half);

		// объединение
		Complex[] freqs = new Complex[n];
		for (int k = 0; k < n/2; k++)
		{
			double kth = -2 * k * Math.PI / n;
			Complex complexExp = new Complex(Math.cos(kth), Math.sin(kth)).multiply(oddFFT[k]);
			freqs[k]        = evenFFT[k].add (complexExp);
			freqs[k + n/2]  = evenFFT[k].add (complexExp).negate();
		}
		return freqs;
	}
}

class AnalyzeData
{
	public static final int CHUNK_SIZE = 512;
	public static final int LOWER_LIMIT = 30;
	public static final int UPPER_LIMIT = 500;
}

class Determinator
{
	public final int[] RANGE = new int[] { 80, 120, 180, 300, AnalyzeData.UPPER_LIMIT+1 };
	private ArrayList<String> hashes = new ArrayList<>();
	private ArrayList<String> freqs = new ArrayList<>();

	// Функция для определения того, в каком диапазоне находится частота
	private int getIndex(int freq)
	{
		int i = 0;
		while (RANGE[i] < freq)
		{
			i++;
		}
		return i;
	}

	public ArrayList<String>[] Determinate(Complex[][] results) throws IOException
	{
		double[] highscores = new double[AnalyzeData.UPPER_LIMIT];
		int[] recordPoints = new int[AnalyzeData.UPPER_LIMIT];

		for (int i = 0; i < results.length; i++)
		{
			for (int freq = AnalyzeData.LOWER_LIMIT; freq < AnalyzeData.UPPER_LIMIT - 1; freq++)
			{
				//Получим силу сигнала
				double mag = Math.log(results[i][freq].abs() + 1);

				//Выясним, в каком мы диапазоне:
				int index = getIndex(freq);

				//Сохраним самое высокое значение силы сигнала и соответствующую частоту
				if (mag > highscores[index])
				{
					highscores[index] = mag;
					recordPoints[index] = freq;
				}
			}
			//Cформируем хэш-тег
			long h = hash(recordPoints[0], recordPoints[1], recordPoints[2], recordPoints[3], recordPoints[4]);
			freqs.add(recordPoints[0] + " " + recordPoints[1] + " " + recordPoints[2] + " " + recordPoints[3] + " " + recordPoints[4]);
			if (h == 0)
			{
				hashes.add(i, "00000000000");
			}
			else
			{
				hashes.add(i, String.valueOf(h));
			}
			highscores = new double[AnalyzeData.UPPER_LIMIT];
			recordPoints = new int[AnalyzeData.UPPER_LIMIT];
		}
		return new ArrayList[]{hashes, freqs};
	}
	private static final int FUZ_FACTOR = 2;

	private long hash(long p1, long p2, long p3, long p4, long p5)
	{
		return ((p5 - (p5 % FUZ_FACTOR)) * 100000000
				+ (p4 - (p4 % FUZ_FACTOR)) * 1000000
				+ (p3 - (p3 % FUZ_FACTOR)) * 10000
				+ (p2 - (p2 % FUZ_FACTOR)) * 100
				+ (p1 - (p1 % FUZ_FACTOR)));
	}
}

class HammingDistance
{
	public static long countDistance(long x, long y)
	{
		long num = x^y;
		long count = 0;
		String str = Long.toBinaryString(num);
		for (int i=0; i < str.length(); i++)
		{
			if (str.charAt(i)=='1')
			{
				count++;
			}
		}
		return count;
	}
}

class SubtractionDistance
{
	public static long countDistance(long x, long y)
	{
		return (long) Math.sqrt(Math.abs(x - y));
	}
}

class MyFrame extends JFrame
{
	MyFrame()
	{
		this.setBounds(0, 0,500,500);
		this.setSize(500,500);
		this.setLayout(null);
		this.setVisible(true);
		this.setLocationRelativeTo(null);
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
}

