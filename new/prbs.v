`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.02.2022 15:53:02
// Design Name: 
// Module Name: prbs
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module prbs(
    input clk, reset,
    
//    output reg [10 : 0] PRBS_OUT, 
//    output wire [10 : 0] PRBS_OUT_BAR
    output reg [2047 : 0] prbs
    );
    localparam IDLE = 1'b0, WORKING = 1'b1;
    
    reg [10 : 0] PRBS_OUT;
    wire [10 : 0] PRBS_OUT_BAR;
    assign PRBS_OUT_BAR = ~PRBS_OUT;
    reg [10 : 0] cnt;
    reg state;
    
    always @ (posedge clk)
    begin
        if(reset == 1) begin
            PRBS_OUT <= 10'b0;
            cnt <= 10'b0;
            state <= WORKING;
        end
        
        else begin 
            if(state == WORKING) begin
                PRBS_OUT[9 : 0] <= PRBS_OUT[10 : 1];
                PRBS_OUT[10] <= ~(PRBS_OUT[0] ^ PRBS_OUT[2]);
                cnt <= cnt + 1;
                prbs[cnt] <= PRBS_OUT[10]; // Assuming prbs11 of matlab is output of 0th bit
                
                if(cnt == 2046) //Finished writing PRBS output to RAM
                    state <= IDLE;
            end
            
        end
    end
    
endmodule
