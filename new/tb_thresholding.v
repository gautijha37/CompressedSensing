`timescale 1ns / 1ps

module tb_thresholding;
//Inputs 
reg clk, reset, enable;
reg signed [11 : 0] sig;

//Outputs
wire signed [11:0] y;
wire done;
wire [11 : 0] id;

////Inputs and outputs of RAM
//wire [9:0] addrb;
//wire signed [31:0] doutb, in_b; wire we_b;

parameter N = 2048;// No. of data samples

////RAM for reading a polynomial values
//poly_dual_file_32_dec  #(.MEM_FILE("input.mif")) read_a(
//            clk,
//            addra,    
//            douta,    
//            addrb,    
//            doutb
//        );

////RAM for writing a0 polynomial values
//ram_32bit_1024_dual store_a0
//(
//    .clk(clk),
//    .in_a(in_0),
//    .in_b(in_b),
//    .addr_a(address_0),
//    .addr_b(addr_b),
//    .we_a(we_0), 
//    .we_b(we_b)
//);

////RAM for writing a1 polynomial values
//ram_32bit_1024_dual store_a1
//(
//    .clk(clk),
//    .in_a(in_1),
//    .in_b(in_b),
//    .addr_a(address_1),
//    .addr_b(addr_b),
//    .we_a(we_1), 
//    .we_b(we_b)
//);

//Thresholding Module Instantiation
thresholding t1(.clk(clk), .reset(reset), .enable(enable),
                .sig(sig),
                
                .y(y),
                .done(done),
                .id(id));              

//Read ecg digital words into ram 
parameter OUT_FILEA0 = "sig_vals.mif";
reg [11 : 0] outa0 [0 : N - 1];

integer fid_out0, j0;
initial begin
    fid_out0 = $fopen (OUT_FILEA0, "r");
    j0 = 0;
    while(! $feof(fid_out0)) begin
        $fscanf(fid_out0, "%d\n", outa0[j0]);
        j0 = j0 + 1;
    end    
    $fclose(fid_out0);
end

//Stimulus
integer addr;
reg test_state;

//Store outputs y into text file
parameter OUT_FILEA1 = "output.txt";
integer fid_out, j;

initial begin
    fid_out = $fopen (OUT_FILEA1, "w");
end

always @ (posedge clk) begin
    if(done == 1) begin
        $fdisplay(fid_out, y);
    end
end

initial
begin
    clk = 1'b0;
    addr = 0;
//    test_state = 0;
    repeat(4000)
    #5 clk = ~clk;
    $fclose(fid_out);
    $stop;
end

initial
begin
    #15 reset <= 1'b1;
    #10 reset <= 1'b0; 
    #10 enable <= 1'b1; 
    #10 enable <= 1'b0;
end

//Sending ecg values from text file to module under test using sig variable
always @ (posedge clk)
begin
    if(reset) begin
        addr <= 0;
//        sig <= outa0[addr];
        test_state <= 1'b0;
    end
    
    else begin
        case(test_state)
        0 : begin
            if(enable) begin
                test_state <= 1;
                sig <= outa0[addr];
                addr <= addr + 1;
            end
        end
        
        1 : begin
            sig <= outa0[addr];
            addr <= addr + 1;
        end
        
        endcase
    end
end

//Comparing output with text file out of C code
//integer it = 0, f0 = 1, f1 = 1;
//always @ (posedge clk) begin
//    $timeformat(-9, 3, " ns");
//    if(done == 1) begin
//        for(it = 0; it < N; it = it + 1) begin
//            if(outa0[it] != store_a0.ram[it]) begin
//                $display("[T = %0t], a0 wrong at i = %0d \n", $realtime, it);
//                f0 = 0;
//            end
            
//            if(outa1[it] != store_a1.ram[it]) begin
//                $display("[T = %0t], a1 wrong at i = %0d \n", $realtime, it);
//                f1 = 0;
//            end
//        end
        
//        if(f0 == 1) begin
//            $display("[T = %0t], a0 correct \n", $realtime);
//        end
        
//        if(f1 == 1) begin
//            $display("[T = %0t], a1 correct \n", $realtime);
//        end
//    end
//end

endmodule
