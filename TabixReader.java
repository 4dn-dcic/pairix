import net.sf.samtools.util.BlockCompressedInputStream;

import java.io.*;
import java.nio.*;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.lang.StringBuffer;

public class TabixReader
{
	private String mFn;
	private BlockCompressedInputStream mFp;

	private int mPreset;
	private int mSc;
	private int mBc;
	private int mEc;
	private int mMeta;
	private int mSkip;
	private String[] mSeq;

	private HashMap<String, Integer> mChr2tid;

	private static int MAX_BIN = 37450;
	private static int TAD_MIN_CHUNK_GAP = 32768;
	private static int TAD_LIDX_SHIFT = 14;

	private class TPair64 implements Comparable<TPair64> {
		long u, v;
		public TPair64(final long _u, final long _v) {
			u = _u; v = _v;
		}
		public TPair64(final TPair64 p) {
			u = p.u; v = p.v;
		}
		public int compareTo(final TPair64 p) {
			return u == p.u? 0 : ((u < p.u) ^ (u < 0) ^ (p.u < 0))? -1 : 1; // unsigned 64-bit comparison
		}
	};

	private class TIndex {
		HashMap<Integer, TPair64[]> b;
		long[] l;
	};
	private TIndex[] mIndex;

	private class TIntv {
		int tid, beg, end;
	};

	private static boolean less64(final long u, final long v) { // unsigned 64-bit comparison
		return (u < v) ^ (u < 0) ^ (v < 0);
	}

	/**
	 * The constructor
	 *
	 * @param fn File name of the data file
	 */
	public TabixReader(final String fn) throws IOException {
		mFn = fn;
		mFp = new BlockCompressedInputStream(new File(fn));
		readIndex();
	}

	private static int reg2bins(final int beg, final int _end, final int[] list) {
		int i = 0, k, end = _end;
		--end;
		list[i++] = 0;
		for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
		for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
		for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
		for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
		for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
		return i;
	}

	public static int readInt(final InputStream is) throws IOException {
		byte[] buf = new byte[4];
		is.read(buf);
		return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getInt();
	}

	public static long readLong(final InputStream is) throws IOException {
		byte[] buf = new byte[8];
		is.read(buf);
		return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getLong();
	}

	public static String readLine(final InputStream is) throws IOException {
		StringBuffer buf = new StringBuffer();
		int c;
		while ((c = is.read()) >= 0 && c != '\n')
			buf.append((char)c);
		if (c < 0) return null;
		return buf.toString();
	}

	/**
	 * Read the Tabix index from a file
	 *
	 * @param fp File pointer
	 */
	public void readIndex(final File fp) throws IOException {
		if (fp == null) return;
		BlockCompressedInputStream is = new BlockCompressedInputStream(fp);
		byte[] buf = new byte[4];

		is.read(buf, 0, 4); // read "TBI\1"
		mSeq = new String[readInt(is)]; // # sequences
		mChr2tid = new HashMap<String, Integer>();
		mPreset = readInt(is);
		mSc = readInt(is);
		mBc = readInt(is);
		mEc = readInt(is);
		mMeta = readInt(is);
		mSkip = readInt(is);
		// read sequence dictionary
		int i, j, k, l = readInt(is);
		buf = new byte[l];
		is.read(buf);
		for (i = j = k = 0; i < buf.length; ++i) {
			if (buf[i] == 0) {
				byte[] b = new byte[i - j];
				System.arraycopy(buf, j, b, 0, b.length);
				String s = new String(b);
				mChr2tid.put(s, k);
				mSeq[k++] = s;
				j = i + 1;
			}
		}
		// read the index
		mIndex = new TIndex[mSeq.length];
		for (i = 0; i < mSeq.length; ++i) {
			// the binning index
			int n_bin = readInt(is);
			mIndex[i] = new TIndex();
			mIndex[i].b = new HashMap<Integer, TPair64[]>();
			for (j = 0; j < n_bin; ++j) {
				int bin = readInt(is);
				TPair64[] chunks = new TPair64[readInt(is)];
				for (k = 0; k < chunks.length; ++k) {
					long u = readLong(is);
					long v = readLong(is);
					chunks[k] = new TPair64(u, v); // in C, this is inefficient
				}
				mIndex[i].b.put(bin, chunks);
			}
			// the linear index
			mIndex[i].l = new long[readInt(is)];
			for (k = 0; k < mIndex[i].l.length; ++k)
				mIndex[i].l[k] = readLong(is);
		}
		// close
		is.close();
	}

	/**
	 * Read the Tabix index from the default file.
	 */
	public void readIndex() throws IOException {
		readIndex(new File(mFn + ".tbi"));
	}

	/**
	 * Read one line from the data file.
	 */
	public String readLine() throws IOException {
		return readLine(mFp);
	}

	private int chr2tid(final String chr) {
		if (mChr2tid.containsKey(chr)) return mChr2tid.get(chr);
		else return -1;
	}

	/**
	 * Parse a region in the format of "chr1", "chr1:100" or "chr1:100-1,000"
	 *
	 * @param reg Region string
	 * @return An array where the three elements are sequence_id,
	 *         region_begin and region_end. On failure, sequence_id==-1.
	 */
	public int[] parseReg(final String reg) { // FIXME: NOT working when the sequence name contains : or -.
		String chr;
		int colon, hyphen;
		int[] ret = new int[3];
		colon = reg.indexOf(':'); hyphen = reg.indexOf('-');
		chr = colon >= 0? reg.substring(0, colon) : reg;
		ret[1] = colon >= 0? Integer.parseInt(reg.substring(colon+1, hyphen)) - 1 : 0;
		ret[2] = hyphen >= 0? Integer.parseInt(reg.substring(hyphen+1)) : 0x7fffffff;
		ret[0] = chr2tid(chr);
		return ret;
	}

	private TIntv getIntv(final String s) {
		TIntv intv = new TIntv();
		int col = 0, end = 0, beg = 0;
		while ((end = s.indexOf('\t', beg)) >= 0) {
			++col;
			if (col == mSc) {
				intv.tid = chr2tid(s.substring(beg, end));
			} else if (col == mBc) {
				intv.beg = intv.end = Integer.parseInt(s.substring(beg, end));
				if ((mPreset&0x10000) != 0) ++intv.end;
				else --intv.beg;
			} else {
				if ((mPreset&0xffff) == 0) { // generic
					if (col == mEc)
						intv.end = Integer.parseInt(s.substring(beg, end));
				} else if ((mPreset&0xffff) == 1) { // SAM
				} else if ((mPreset&0xffff) == 2) { // VCF
				}
			}
			beg = end + 1;
		}
		return intv;
	}

	public class Iterator {
		private int i, n_seeks;
		private int tid, beg, end;
		private TPair64[] off;
		private long curr_off;
		private boolean iseof;

		public Iterator(final int _tid, final int _beg, final int _end, final TPair64[] _off) {
			i = -1; n_seeks = 0; curr_off = 0; iseof = false;
			off = _off; tid = _tid; beg = _beg; end = _end;
		}

		public String next() throws IOException {
			if (iseof) return null;
			for (;;) {
				if (curr_off == 0 || !less64(curr_off, off[i].v)) { // then jump to the next chunk
					if (i == off.length - 1) break; // no more chunks
					if (i >= 0) assert(curr_off == off[i].v); // otherwise bug
					if (i < 0 || off[i].v != off[i+1].u) { // not adjacent chunks; then seek
						mFp.seek(off[i+1].u);
						curr_off = mFp.getFilePointer();
						++n_seeks;
					}
					++i;
				}
				String s;
				if ((s = readLine(mFp)) != null) {
					TIntv intv;
					char[] str = s.toCharArray();
					curr_off = mFp.getFilePointer();
					if (str.length == 0 || str[0] == mMeta) continue;
					intv = getIntv(s);
					if (intv.tid != tid || intv.beg >= end) break; // no need to proceed
					else if (intv.end > beg && intv.beg < end) return s; // overlap; return
				} else break; // end of file
			}
			iseof = true;
			return null;
		}
	};

	public Iterator query(final int tid, final int beg, final int end) {
		TPair64[] off, chunks;
		long min_off;
		TIndex idx = mIndex[tid];
		int[] bins = new int[MAX_BIN];
		int i, l, n_off, n_bins = reg2bins(beg, end, bins);
		min_off = (beg>>TAD_LIDX_SHIFT >= idx.l.length)? 0 : idx.l[beg>>TAD_LIDX_SHIFT];
		for (i = n_off = 0; i < n_bins; ++i) {
			if ((chunks = idx.b.get(bins[i])) != null)
				n_off += chunks.length;
		}
		if (n_off == 0) return null;
		off = new TPair64[n_off];
		for (i = n_off = 0; i < n_bins; ++i)
			if ((chunks = idx.b.get(bins[i])) != null)
				for (int j = 0; j < chunks.length; ++j)
					if (less64(min_off, chunks[j].v))
						off[n_off++] = new TPair64(chunks[j]);
		Arrays.sort(off, 0, n_off);
		// resolve completely contained adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i) {
			if (less64(off[l].v, off[i].v)) {
				++l;
				off[l].u = off[i].u; off[l].v = off[i].v;
			}
		}
		n_off = l + 1;
		// resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
		for (i = 1; i < n_off; ++i)
			if (!less64(off[i-1].v, off[i].u)) off[i-1].v = off[i].u;
		// merge adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i) {
			if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
			else {
				++l;
				off[l].u = off[i].u;
				off[l].v = off[i].v;
			}
		}
		n_off = l + 1;
		// return
		TPair64[] ret = new TPair64[n_off];
		for (i = 0; i < n_off; ++i) ret[i] = new TPair64(off[i].u, off[i].v); // in C, this is inefficient
		return new TabixReader.Iterator(tid, beg, end, ret);
	}
	
	public Iterator query(final String reg) {
		int[] x = parseReg(reg);
		return query(x[0], x[1], x[2]);
	}

	public static void main(String[] args) {
		if (args.length < 1) {
			System.out.println("Usage: java -cp .:sam.jar TabixReader <in.gz> [region]");
			System.exit(1);
		}
		try {
			TabixReader tr = new TabixReader(args[0]);
			String s;
			if (args.length == 1) { // no region is specified; print the whole file
				while ((s = tr.readLine()) != null)
					System.out.println(s);
			} else { // a region is specified; random access
				TabixReader.Iterator iter = tr.query(args[1]); // get the iterator
				while ((s = iter.next()) != null)
					System.out.println(s);
			}
		} catch (IOException e) {
		}
	}
}
