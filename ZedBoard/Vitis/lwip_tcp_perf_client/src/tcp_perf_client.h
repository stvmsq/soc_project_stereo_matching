/*
 * Copyright (C) 2018 - 2019 Xilinx, Inc.
 * Copyright (C) 2022 - 2024 Advanced Micro Devices, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 * OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

#ifndef __TCP_PERF_CLIENT_H_
#define __TCP_PERF_CLIENT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "lwipopts.h"
#include "lwip/ip_addr.h"
#include "lwip/err.h"
#include "lwip/tcp.h"
#include "lwip/inet.h"
#include "xil_printf.h"
#include "platform.h"

#include "frame_buffer.h"

/*enum {
	KCONV_UNIT,
	KCONV_KILO,
	KCONV_MEGA,
	KCONV_GIGA,
};

const char kLabel[] =
{
	' ',
	'K',
	'M',
	'G'
};

enum measure_t {
	BYTES,
	SPEED
};

enum report_type {
	INTER_REPORT,
	TCP_DONE_CLIENT,
	TCP_ABORTED_REMOTE
};

struct interim_report {
	u64_t start_time;
	u64_t last_report_time;
	u32_t total_bytes;
	u32_t report_interval_time;
};

struct perf_stats {
	u8_t client_id;
	u64_t start_time;
	u64_t end_time;
	u64_t total_bytes;
	struct interim_report i_report;
};*/

typedef enum {
    WAIT,      // 0
	WAIT_FOR_FRAME,      // 1
    RECEIVE_1,   // 2
	RECEIVE_2,	// 3
	SEND_IMG,	// 4
	NOTIFY_TO_CLOSE, // 5
    CLOSE,       // 6
} ConnectionStatus_t;

void start_tcp_client(SharedMemory_t *frame_buffer);
err_t send_depth_img(void);
err_t request_new_frame(void);
ConnectionStatus_t get_connection_status(void);
void tcp_connection_close(void);
/* seconds between periodic bandwidth reports */
//#define INTERIM_REPORT_INTERVAL 5

/* Client port to connect */
#define TCP_CONN_PORT 5001

/* time in seconds to transmit packets */
#define TCP_TIME_INTERVAL 300

#if LWIP_IPV6==1
/* Server to connect with */
#define TCP_SERVER_IPV6_ADDRESS "fe80::6600:6aff:fe71:fde6"
#else
/* Server to connect with */
#define TCP_SERVER_IP_ADDRESS "192.168.1.20"
#endif

// #define TCP_SEND_BUFSIZE (5*TCP_MSS)

#ifdef __cplusplus
}
#endif

#endif /* __TCP_PERF_CLIENT_H_ */
