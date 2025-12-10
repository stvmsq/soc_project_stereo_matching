# 0 "/home/xilinx/Projects/stere_matching_software/vitis_software/platform/zynq_fsbl/zynq_fsbl_bsp/lop-config.dts"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/home/xilinx/Projects/stere_matching_software/vitis_software/platform/zynq_fsbl/zynq_fsbl_bsp/lop-config.dts"

/dts-v1/;
/ {
        compatible = "system-device-tree-v1,lop";
        lops {
                lop_0 {
                        compatible = "system-device-tree-v1,lop,load";
                        load = "assists/baremetal_validate_comp_xlnx.py";
                };

                lop_1 {
                    compatible = "system-device-tree-v1,lop,assist-v1";
                    node = "/";
                    outdir = "/home/xilinx/Projects/stere_matching_software/vitis_software/platform/zynq_fsbl/zynq_fsbl_bsp";
                    id = "module,baremetal_validate_comp_xlnx";
                    options = "ps7_cortexa9_0 /tools/Xilinx/2025.1/Vitis/data/embeddedsw/lib/sw_services/xilffs_v5_4/src /home/xilinx/Projects/stere_matching_software/vitis_software/_ide/.wsdata/.repo.yaml";
                };

                lop_2 {
                    compatible = "system-device-tree-v1,lop,assist-v1";
                    node = "/";
                    outdir = "/home/xilinx/Projects/stere_matching_software/vitis_software/platform/zynq_fsbl/zynq_fsbl_bsp";
                    id = "module,baremetal_validate_comp_xlnx";
                    options = "ps7_cortexa9_0 /tools/Xilinx/2025.1/Vitis/data/embeddedsw/lib/sw_services/xilrsa_v1_8/src /home/xilinx/Projects/stere_matching_software/vitis_software/_ide/.wsdata/.repo.yaml";
                };

        };
    };
