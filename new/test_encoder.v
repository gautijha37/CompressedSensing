`timescale 1ns / 1ps

module test_encoder;

//thresholding module inputs
reg clk, reset, enable;
reg signed [11 : 0] sig;

//thresholding module outputs
wire signed [11:0] y;
wire done;
wire [11 : 0] id;

parameter N = 2048;// No. of data samples

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

//Thresholding Module Instantiation
thresholding t1(.clk(clk), .reset(reset), .enable(enable),
                .sig(sig),
                
                .y(y),
                .done(done),
                .id(id));    

//Stimulus
integer addr;
reg test_state;
        
//Stimulus
initial
begin
    clk = 1'b0;
    addr = 0;
    repeat(5000)
    #5 clk = ~clk;
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

////Instantiating prbs module
//wire [2047 : 0] prbs;
//reg prbs_reset;
//prbs p1(.clk(clk), .reset(prbs_reset), 
//        .prbs(prbs));

//Read ecg digital words into ram 
parameter prbsfile = "prbsvals.mif";
reg [2047 : 0] prbs_in;

integer fid_prbs, j3;
initial begin
    fid_prbs = $fopen (prbsfile, "r");
    j3 = 0;
    while(! $feof(fid_prbs)) begin
        $fscanf(fid_prbs, "%d\n", prbs_in[j3]);
        j3 = j3 + 1;
    end    
    $fclose(fid_prbs);
end

wire encoder_done;

//Encoder instantiation
encoder e1(.clk(clk), .reset(reset), .doneThreshold(done),
           .threshold_sig(y),
           .prbs_in(prbs_in),
           .id(id),
           
           .done(encoder_done));

//Writing final output to text file
parameter OUT_FILE = "final_out.txt";
integer fid_out, i;

initial begin
    fid_out = $fopen (OUT_FILE, "w");
end

always @ (posedge clk) begin
    if(id == 2047) begin
        for(i = 0; i < 512; i = i + 1) 
            $fdisplay(fid_out, e1.y[i]);
    end
end
endmodule
