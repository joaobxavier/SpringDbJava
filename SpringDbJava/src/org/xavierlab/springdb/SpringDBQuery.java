package org.xavierlab.springdb;

import java.io.BufferedReader;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;

public class SpringDBQuery {
	private String url;
	private Connection connection;

	/**
	 * Connects to the database
	 * 
	 * @throws SQLException
	 */
	public SpringDBQuery(String dbConfigFile) throws SQLException {
		// use this to use a local version of SpringDb
		// this.useLocalDb();
		// use this to use the cloud version
		this.connectToDb(dbConfigFile);
	}

	/**
	 * Make the connection to the database
	 * 
	 * @throws SQLException
	 */
	private void connectToDb(String dbConfigFile) throws SQLException {
		// connect to the database
		buildUrlFromConfigFile(dbConfigFile);
		connection = DriverManager.getConnection(url);
		System.out.println(this.toString());
	}

	private void buildUrlFromConfigFile(String dbConfigFile) {
		try {
			FileReader fr = new FileReader(dbConfigFile);
			BufferedReader textReader = new BufferedReader(fr);
			String textData = textReader.readLine();
			textReader.close();
			String[] tokens = textData.split(" ");
			// get the address (first token)
			String address = (tokens[0].split("'"))[1];
			// get the database (2nd token)
			String database = (tokens[1].split("'"))[1];
			// get the user (3rd token)
			String user = (tokens[2].split("'"))[1];
			// get the password (4th token)
			String password = (tokens[3].split("'"))[1];
			// concatenate
			url = "jdbc:postgresql://" + address + "/" + database + "?"
					+ "user=" + user + "&password=" + password + "&ssl=true"
					+ "&sslfactory=org.postgresql.ssl.NonValidatingFactory";
			System.out.println(url);

		} catch (Exception e) {
			System.out
					.println("must provide config file with the follwoing info:");
			System.out
					.println("host='#####' dbname='#####' user='#####' password='#####'");
			e.printStackTrace();
		}
	}

	/**
	 * Runs an SQL query, returns the result as a matrix of strings
	 * 
	 * @param query
	 * @return
	 * @throws SQLException
	 */
	public String[][] returnQuery(String query) throws SQLException {

		Statement st = connection.createStatement(
				ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_UPDATABLE);

		// get genomes
		ResultSet rs = st.executeQuery(query);
		ResultSetMetaData rsmd = rs.getMetaData();
		int numberOfColums = rsmd.getColumnCount();
		// count number of rows
		int count = 0;
		while (rs.next()) {
			++count;
			// Get data from the current row and use it
		}
		if (count == 0) {
			System.out.println("No records found");
		}
		// return to begining
		rs.beforeFirst();
		// initialize table
		String[][] names = new String[count][numberOfColums];
		for (int i = 0; i < count; i++) {
			rs.next();
			for (int j = 0; j < numberOfColums; j++) {
				String name = rs.getString(j + 1);
				names[i][j] = name;
			}
		}
		return names;
	}

	public float[] getPhenotypeColumn(String phenotype) throws SQLException {
		String[][] r = this
				.returnQuery("SELECT "
						+ phenotype
						+ " from"
						+ " phenotype where genome_id in"
						+ " (SELECT DISTINCT genome_id FROM orf WHERE orth_orf_id IS NOT NULL)"
						+ " and genome_id in (SELECT genome_id FROM phenotype) ORDER BY genome_id");
		return convertColumnToArrayOfNumbers(r, 0);
	}

	public String getOrth_Orf_Id_OfGene(String geneName) throws SQLException {
		String[][] r = returnQuery("select distinct(orth_orf_id) from orf where gene_name like '%"
				+ geneName + "%'");
		return r[0][0];
	}

	public String[] getSequencesOfGene(String gene) throws SQLException {
		String[][] r = this
				.returnQuery("SELECT seq from"
						+ " orf where genome_id in"
						+ " (SELECT DISTINCT genome_id FROM orf WHERE orth_orf_id IS NOT NULL)"
						+ " and orth_orf_id="
						+ getOrth_Orf_Id_OfGene(gene)
						+ " and genome_id in (SELECT genome_id FROM phenotype) ORDER BY genome_id");
		String[] r2 = new String[r.length];
		for (int i = 0; i < r2.length; i++) {
			r2[i] = r[i][0];
		}
		return r2;

	}

	/**
	 * Same as getSequencesOfGene but can request for a given genome
	 * 
	 * @param gene
	 * @param genome_id
	 * @return
	 * @throws SQLException
	 */
	public String getSequencesOfGene(String gene, int genome_id)
			throws SQLException {
		String[][] r = this.returnQuery("SELECT seq from"
				+ " orf where orth_orf_id=" + getOrth_Orf_Id_OfGene(gene)
				+ " and genome_id = " + genome_id);
		return r[0][0];
	}

	public String[][] getSequencesOfGenes(String[] genes) throws SQLException {
		// construct sql query
		String q = "SELECT B.gene_name, B.orth_orf_id, A.genome_id, A.seq from"
				+ " orf A RIGHT JOIN (SELECT distinct orth_orf_id, gene_name from orf where ";
		for (int i = 0; i < genes.length; i++) {
			q += "gene_name like '%" + genes[i] + "%'";
			if (i < (genes.length - 1)) {
				q += " or ";
			}	
		}
		q += ") B on A.orth_orf_id = B.orth_orf_id and genome_id in (SELECT genome_id FROM phenotype) ORDER BY genome_id";
		//
		String[][] r = this.returnQuery(q);
		// get the orth_orf_ids for the genes (in string format)
		String orth_orf_ids[] = new String[genes.length];
		for (int i = 0; i < genes.length; i++) {
			for (int j = 0; j < r.length; j++) {
				if (r[j][0].contains(genes[i])){
					orth_orf_ids[i] = r[j][1];
					break;
				}		
			}
		}
		// get the genomeIds
		int[] genome_ids_in_query = new int[r.length];
		HashSet<String> tmp1 = new HashSet<String>();
		for (int i = 0; i < r.length; i++) {
			tmp1.add(r[i][2]);
			genome_ids_in_query[i] = Integer.parseInt(r[i][2]);
		}
		String[] tmp2 = new String[tmp1.size()];
		tmp2 = tmp1.toArray(tmp2);
		int[] genome_ids = new int[tmp2.length];
		for (int i = 0; i < tmp2.length; i++) {
			genome_ids[i] = Integer.parseInt(tmp2[i]);
		}
		Arrays.sort(genome_ids);
		// prepare the output: a matrix of String objects where each entry is a sequence
		// the matrix is ordered so that rows represent genomes (sorted in ascending order of genoem_id)
		// and columns are ordered following the input list
		String[][] sequences = new String[genome_ids.length][orth_orf_ids.length];
		for (int i = 0; i < r.length; i++) {
			int row = 0;
			for (int j = 0; j < genome_ids.length; j++) {
				if (genome_ids_in_query[i] == genome_ids[j]) {
					row = j;
					break;
				}
			}
			int column = 0;
			for (int j = 0; j < orth_orf_ids.length; j++) {
				if (r[i][1].equals(orth_orf_ids[j])) {
					column = j;
					break;
				}
			}
			sequences[row][column] = r[i][3];	
		}
		return sequences;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return "Connected to " + url;
	}

	public static float[] convertColumnToArrayOfNumbers(String[][] s,
			int columnNumber) {
		float[] n = new float[s.length];
		for (int i = 0; i < s.length; i++) {
			n[i] = Float.parseFloat(s[i][columnNumber]);
		}
		return n;
	}

	/**
	 * Main method can be used for running the class as a command or for test
	 * purpuses
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		try {
			SpringDBQuery dB = new SpringDBQuery(args[0]);
			// sprindDbQuery.useLocalDb();
			String[][] result;
			result = dB
					.returnQuery("SELECT genome_id, strain_name from"
							+ " genome where genome_id in"
							+ " (SELECT DISTINCT genome_id FROM orf WHERE orth_orf_id IS NOT NULL)"
							+ " and genome_id in (SELECT genome_id FROM phenotype) ORDER BY genome_id");
			for (int i = 0; i < result.length; i++) {
				for (int j = 0; j < result[0].length; j++) {
					String s = result[i][j];
					System.out.print(s + "\t");
				}
				System.out.println();
			}
			float[] r = dB.getPhenotypeColumn("biofilm");
			for (int i = 0; i < r.length; i++) {
				System.out.println(r[i]);
			}
			System.out.println(dB.getOrth_Orf_Id_OfGene("wspF"));
			String[] seqs = dB.getSequencesOfGene("wspR");
			for (int i = 0; i < seqs.length; i++) {
				System.out.println(seqs[i]);
			}
			System.out.println("\n\nTest getting a single gene from PA14");
			System.out.println(dB.getSequencesOfGene("wspF", 11));

			//
			String[] genes = new String[2];
			genes[0] = "wspF";
			genes[1] = "wspE";
			String[][] seqMatrix = dB.getSequencesOfGenes(genes);
			for (int i = 0; i < seqMatrix.length; i++) {
				for (int j = 0; j < seqMatrix[i].length; j++) {
					System.out.println(seqMatrix[i][j]);					
				}
				System.out.println("\n\n");
			}

		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
