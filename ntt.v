`include "fmt.v"

/*Main compute module, includes everything for quick memory access*/

module compute(
	input clk, 
	input rst, 
	input [`WORD-1:0] cycles,
	output reg done
	);

	integer ind;
	
	// general
    reg [`DWORD-1:0] conv_size;
    reg [`WORD-1:0] rc1, rc2;

    reg [`DWORD-1:0] x [`MAXOUT-1:0];
    reg [`DWORD-1:0] y[`MAXOUT-1:0];

    reg [`DWORD-1:0] roots [`MAXSIZE:0];
    reg [`DWORD-1:0] inv_roots [`MAXSIZE:0];
    reg [`DWORD-1:0] inv_conv [`MAXSIZE:0];

    reg [`WORD-1:0] rev [`MAXSIZE:0];

    reg [`DWORD-1:0] out [`MAXOUT-1:0];

    // comp enable regs
	reg convp_comp, prerev_comp, reverse_comp, NTT_comp, dot_comp, inv_comp, norm_comp;
    // status flags
    reg conv_f, roots_f, NTT_f, prerev_f, reverse_f, dot_f, inv_f, norm_f;

    // repeat process
    reg revcnt, nttcnt;

	// compute convolution size
	wire [`DWORD-1:0] raw_size;
	assign raw_size = rc1+rc2-1;
	reg [`DWORD-1:0] log_size;
	always @(posedge clk) begin
		if(convp_comp) begin
			conv_size <= 2**($clog2(raw_size));
			log_size <= $clog2(raw_size);
			convp_comp <= 0;
			conv_f <= 1;
		end
	end

	// precompute bit reversal locations
	reg [`DWORD-1:0] prerev_ind, cur;
	reg [`DWORD-1:0] j;
	wire tr1, tr2;
	assign tr1 = cur[j];
	assign tr2 = cur[log_size-j-1];
	always @(posedge clk) begin
		if(prerev_comp) begin
			if(loop==0) begin
				// $display("%d", prerev_ind);
				if(prerev_ind<conv_size) begin
					// $display("in here");
					j <= 0;
					loop <= 1;
					cur <= prerev_ind;
				end
				else begin
					prerev_ind <= 0;
					cur <= 0;
					prerev_comp <= 0;
					prerev_f <= 1;
				end
			end
			else begin
				if(j<log_size/2) begin
					cur[j] <= tr2;
					cur[log_size-j-1] <= tr1;
					j<=j+1;
				end
				else begin
					rev[prerev_ind] <= cur;
					// $display("HERE %d", cur);
					loop <= 0;
					prerev_ind <= prerev_ind + 1;
				end
			end
			
		end
	end

	// reverse bits
	wire [`DWORD-1:0] half;
	assign half = conv_size/2;
	reg [`DWORD-1:0] rev_ind;
	reg dummy;

	wire [`WORD-1:0] tx1, tx2, ty1, ty2, to1, to2; // for swapping
	assign tx1 = x[rev_ind];
	assign tx2 = x[rev[rev_ind]];
	assign ty1 = y[rev_ind]; 
	assign ty2 = y[rev[rev_ind]];

	assign to1 = out[rev_ind];
	assign to2 = out[rev[rev_ind]];

	always @(posedge clk) begin
		if(reverse_comp) begin
			if(rev_ind < conv_size) begin
				if(rev_ind<rev[rev_ind]) begin
					x[rev_ind] <= tx2;
					x[rev[rev_ind]] <= tx1;
					y[rev_ind] <= ty2;
					y[rev[rev_ind]] <= ty1;
					out[rev_ind] <= to2;
					out[rev[rev_ind]] <= to1;
				end
				rev_ind <= rev_ind + 1;
			end
			else begin
				rev_ind <= 0;
				reverse_comp <= 0;
				reverse_f <= 1;
			end
		end
	end


	// NTT 
	reg [`DWORD-1:0] len;
	reg [`DWORD-1:0] omega, w;
	reg [1:0] loop;
	reg [`DWORD-1:0] i, k;
	reg [`DWORD-1:0] u, v;
	reg [`DWORD-1:0] half_ntt;
	always @(posedge clk) begin
		if(NTT_comp) begin
			$display("doing NTT");
			if(loop==0) begin
				if((len)<=conv_size) begin
					if(inv_comp) omega <= inv_roots[len];
					else omega <= roots[len];
					if(!inv_comp) $display("omega: %d", roots[len]);
					// len <= len<<1;
					loop <= 1;
					i <= 0;
					k <= 0;
				end
				else begin
					NTT_comp <= 0;
					NTT_f <= 1;
					len <= 2;
				end
			end
			else if(loop==1) begin
				// $display("%d, %d", i, conv_size);
				if(i<conv_size) begin
					// $display("SKIBIDI <<<<<<<<");
					half_ntt <= len/2;
					w <= 1;
					k <= 0;
					loop <= 2;
					// i <= i + len;
				end
				else begin
					loop<=0;
					len <= len<<1;
				end
			end
			else if(loop==2) begin
				if(k<half_ntt) begin
					// $display("SIGMA --------------");
					x[i+k] <= (x[i+k] + ((x[i+k+half_ntt]*w)%`MODP))%`MODP;
					x[i+k+half_ntt] <= (x[i+k]-((x[i+k+half_ntt]*w)%`MODP)+`MODP)%`MODP;

					y[i+k] <= (y[i+k] + ((y[i+k+half_ntt]*w)%`MODP))%`MODP;
					y[i+k+half_ntt] <= (y[i+k]-((y[i+k+half_ntt]*w)%`MODP)+`MODP)%`MODP;

					out[i+k] <= (out[i+k] + ((out[i+k+half_ntt]*w)%`MODP))%`MODP;
					out[i+k+half_ntt] <= (out[i+k]-((out[i+k+half_ntt]*w)%`MODP)+`MODP)%`MODP;

					k <= k+1;
					w <= (w*omega)%`MODP;
				end
				else begin
					loop<=1;
					i<= i+len;
				end
			end
		end
	end

	// dot product x and y
	reg [`WORD-1:0] dot_ind;
	always @(posedge clk) begin
		if(dot_comp) begin
			if(dot_ind<conv_size) begin
				$display("dot_ind: %d", dot_ind);
				out[dot_ind] <= (x[dot_ind]*y[dot_ind])%`MODP;
				dot_ind <= dot_ind+1;
			end
			else begin
				dot_ind <= 0;
				dot_comp <= 0;
				dot_f <= 1;
			end
		end
	end

	// inverse NTT
	always @(posedge clk) begin
		if(inv_comp) begin
			NTT_comp <= 1;
			if(NTT_f) begin
				inv_comp <= 0;
				inv_f <= 1;
				NTT_f <= 0;
				NTT_comp <= 0;
			end
		end
	end

	// normalize
	reg [`DWORD-1:0] norm_ind;
	always @(posedge clk) begin
		if(norm_comp) begin
			$display("computing skibidi");
			if(norm_ind<conv_size) begin
				$display("computing skibidi: %d", out[norm_ind]);
				out[norm_ind] <= (out[norm_ind]*inv_conv[conv_size])%`MODP;
				norm_ind <= norm_ind + 1;
			end
			else begin
				norm_comp <= 0;
				norm_ind <= 0;
				norm_f <= 1;
			end
		end
	end

	// main control
	reg NTT_cnt;
	always @(posedge clk) begin
		if(rst) begin
			$display("finding convolution size\n");
			convp_comp <= 1;
			NTT_cnt <= 0;
			revcnt <= 0;
		end
		else if(conv_f) begin
			$display("reversing bits\n");
			prerev_comp <= 1;
			conv_f <= 0;
			$display("%16b\n", conv_size);
			$display("%16b\n", log_size);
		end
		if(prerev_f) begin
			$display("precomputing bit reversals\n");
			prerev_f <= 0;
			reverse_comp <= 1;
		end
		else if(reverse_f) begin
			if(revcnt==0) begin
				for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", x[ind]);
				for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", y[ind]);
				$display("doing NTT forward\n");
				revcnt <= 1;
				NTT_comp <= 1;
			end
			else begin
				$display("After 2nd inverse: ");
				for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", out[ind]);
				$display("doing NTT inverse\n");
				revcnt <= 0;
				inv_comp <= 1;
			end
			reverse_f <= 0;
		end
		else if(NTT_f && NTT_cnt==0) begin
			$display("NTT FORWARD RESULT");
			for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", x[ind]);
			for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", y[ind]);
			$display("doing dot product\n");
			dot_comp <= 1;
			NTT_f <= 0;
			NTT_cnt <= 1;
		end
		else if(dot_f) begin
			dot_f <= 0;
			reverse_comp <= 1;
			$display("DOT PRODUCT RESULT:");
			for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", out[ind]);
			$display("reversing bits\n");
		end
		else if(inv_f) begin
			norm_comp <= 1;
			inv_f <= 0;
			for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", out[ind]);
			$display("normalizing: ");			
		end
		else if(norm_f) begin
			$display("%d", inv_conv[conv_size]);
			norm_f <= 0;
			done <= 1;
			$display("res: ");
			for(ind = 0; ind<conv_size; ind = ind+1) $display("%d ", out[ind]);
			$display("clock cycles: %d\n", cycles);
		end
	end


	

	initial begin
		convp_comp <= 0;
		reverse_comp <= 0;
		NTT_comp <= 0;
		loop <= 0;
		len <= 2;
		prerev_ind <= 0;
		cur <= 0;
		rev_ind <= 0;
		dot_ind <= 0;
		dot_comp <= 0;
		inv_comp <= 0;
		norm_comp <= 0;
		norm_ind <= 0;
		revcnt <= 0;
		nttcnt <= 0;
		inv_f <= 0;
		NTT_cnt <= 0;
		k <= 1;
		$readmemh("roots.txt", roots);
		$readmemh("inv_roots.txt", inv_roots);
		$readmemh("inv_n.txt", inv_conv);

		for(ind = 0; ind<`MAXOUT; ind = ind+1) begin
			out[ind] <= 0;
			x[ind] <= 0;
			y[ind] <= 0;
			// rev[ind] <= 0
		end
		ind <= 0;

		// test

		x[0] <= 1;
		x[1] <= 1;
		y[0] <= 0;
		y[1] <= 1;
		rc1 <= 2;
		rc2 <= 2;
	end

endmodule