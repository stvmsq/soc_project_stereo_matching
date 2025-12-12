#include "tcp_perf_client.h"
#include <stdint.h>
#include <xil_printf.h>
#include "frame_buffer.h"
#include "lwip/err.h"
#include "lwip/tcp.h"

static struct tcp_pcb *c_pcb;

ConnectionStatus_t connection_status;
SharedMemory_t *shared_memory;

ConnectionStatus_t get_connection_status(void){
	return connection_status;
}
/** Close a tcp session */
static void tcp_client_close(struct tcp_pcb *pcb)
{
	err_t err;
	xil_printf("INFO: Connection close\r\n");
	if (pcb != NULL) {
		tcp_sent(pcb, NULL);
		tcp_err(pcb, NULL);
		err = tcp_close(pcb);
		if (err != ERR_OK) {
			/* Free memory with abort */
			tcp_abort(pcb);
		}
	}
}

void tcp_connection_close(void){
	tcp_client_close(c_pcb);
}

/** Error callback, tcp session aborted */
static void tcp_client_err(void *arg, err_t err)
{
	if (err == ERR_RST){
		xil_printf("ERROR: Connection reset! Is the Server running with IP: %s and Port: %d ?\n\r", TCP_SERVER_IP_ADDRESS, TCP_CONN_PORT);
	}else{
		xil_printf("ERROR: TCP error: %s (%d)\r\n", lwip_strerr(err), err);
	}
	// LWIP_UNUSED_ARG(err);
	
	tcp_client_close(c_pcb);
	c_pcb = NULL;
	xil_printf("INFO: TCP connection aborted\n\r");
}

static err_t tcp_send_message(const uint8_t *buffer, uint16_t length){
	err_t err;
	if (c_pcb == NULL) {
		return ERR_CONN;
	}
	if (tcp_sndbuf(c_pcb) < length) {
		xil_printf("WARNING: Send buffer full, cannot write!\r\n");
		return ERR_MEM;
	}
	err = tcp_write(c_pcb, buffer, length, TCP_WRITE_FLAG_COPY);
	if (err != ERR_OK) {
		xil_printf("ERROR: TCP client: Error on tcp_write: %d\r\n", err);
		return err;
	}
	err = tcp_output(c_pcb);
	if (err != ERR_OK) {
		xil_printf("ERROR: TCP client: Error on tcp_output: %d\r\n", err);
		return err;
	}
	return ERR_OK;
}

err_t request_new_frame(void){
	if(connection_status != WAIT) {
		xil_printf("WARNING: request new frame during operation: %d\r\n", connection_status);
		return ERR_INPROGRESS;
	}
	if(shared_memory->stereo_frame[next_frame_id(shared_memory)].status == INVALID || shared_memory->stereo_frame[next_frame_id(shared_memory)].status == DONE){
		uint8_t buffer[1];
		buffer[0] = 1;
		if (tcp_send_message(buffer, sizeof(buffer)) == ERR_OK) {
			//ToDO: notify user
			connection_status = WAIT_FOR_FRAME;
		}
	} else {
		xil_printf("WARNING: request new frame, but old frame is not invalid or done.\n\r");
	}
	return ERR_OK;
}

int32_t send_index = -1;
err_t send_depth_img(void){
	if (connection_status != WAIT && connection_status != SEND_IMG) {
		 return ERR_INPROGRESS;
	}
	if (connection_status == WAIT){
		if (send_index == -1 && shared_memory->last_depth_frame.status == READY){
			// send_index and status condition should alway be fuilfilled here
			connection_status = SEND_IMG;
		} else{
			xil_printf("WARING: call send_depth_img in wrong state: %d and status: %d \r\n", send_index, shared_memory->last_depth_frame.status);
			return ERR_VAL;
		}
	}
	if (shared_memory->last_depth_frame.status == READY) {
		if (send_index == -1 && tcp_sndbuf(c_pcb) >= 8) {
			// Buffer for 1B (type) + 1 x uint32 + 2x uint16 = 9 Bytes
			uint8_t header[9];
			header[0] = 3;
			header[1] = (uint8_t)(shared_memory->last_depth_frame.frame_id);
			header[2] = (uint8_t)(shared_memory->last_depth_frame.frame_id >> 8);
			header[3] = (uint8_t)(shared_memory->last_depth_frame.frame_id >> 16);
			header[4] = (uint8_t)(shared_memory->last_depth_frame.frame_id >> 24);

			header[5] = (uint8_t)(shared_memory->last_depth_frame.width);
			header[6] = (uint8_t)(shared_memory->last_depth_frame.width >> 8);

			header[7] = (uint8_t)(shared_memory->last_depth_frame.height);
			header[8] = (uint8_t)(shared_memory->last_depth_frame.height >> 8);

			if (tcp_send_message(header, sizeof(header)) == ERR_OK) {
				send_index = 0;
			}
		} else if (send_index >= 0 && send_index < IMG_HEIGHT * IMG_WIDTH * 4) {
			uint16_t sndbuf = tcp_sndbuf(c_pcb);
			if (sndbuf > 0) {
				// number of bytes we can send now
        		uint16_t chunk = (IMG_HEIGHT * IMG_WIDTH * 4 - send_index < sndbuf) ? IMG_HEIGHT * IMG_WIDTH * 4 - send_index : sndbuf;
				const uint8_t *ptr = (const uint8_t*)shared_memory->last_depth_frame.depth + send_index;
				if (tcp_send_message(ptr, chunk) == ERR_OK) {
					send_index += chunk;
				}
			}
		} else if (send_index >= IMG_HEIGHT * IMG_WIDTH * 4) {
			xil_printf("INFO: send depth img: %d\r\n", shared_memory->last_depth_frame.frame_id);
			shared_memory->last_depth_frame.status = DONE;
			connection_status = WAIT;
			send_index = -1;
		}
	}
	return ERR_OK;
}


/** TCP sent callback, try to send more data */
static err_t tcp_client_sent(void *arg, struct tcp_pcb *tpcb, u16_t len)
{
	return 0;
	// return tcp_send_perf_traffic();
}

static uint32_t rec_index = 1;

static err_t recv_img(void *arg, struct tcp_pcb *tpcb, struct pbuf *p, err_t err) {
	struct pbuf *q;
	// Durch alle pbufs iterieren
    for (q = p; q != NULL; q = q->next) {
        uint8_t *payload = (uint8_t *)q->payload;
		for (uint16_t i = 0; i < q->len;) {
			uint16_t byte_left = q->len - i;
			if (rec_index == 1 && byte_left >= 8) {
				// xil_printf("INFO: copy parameter\n\r");
				// 4B frame_id + 2B width + 2B height
				memcpy(&shared_memory->stereo_frame[next_frame_id(shared_memory)].frame_id, &payload[i], 8);
				rec_index += 8;
				i+=8;
				if (shared_memory->stereo_frame[next_frame_id(shared_memory)].width != IMG_WIDTH || shared_memory->stereo_frame[next_frame_id(shared_memory)].height != IMG_HEIGHT) {
					xil_printf("ERROR: Wrong image size, schould %dx%d but is %dx%d\n\r", IMG_WIDTH, IMG_HEIGHT, (int)shared_memory->stereo_frame[next_frame_id(shared_memory)].width, (int)shared_memory->stereo_frame[next_frame_id(shared_memory)].height);
				}
				if(connection_status == RECEIVE_2){
					// xil_printf("INFO: skip calib\n\r");
					// if no calib data send (type = 2), skip next rec_index
					rec_index += 80;
				}
			}else if (rec_index == 9 && connection_status == RECEIVE_1 && byte_left >= 80) {
				// xil_printf("INFO: copy calib\n\r");
				// 80B calib
				memcpy(&shared_memory->stereo_frame[next_frame_id(shared_memory)].calib, &payload[i], 80);
				rec_index += 80;
				i += 80;
			}else if (rec_index >= 89 && rec_index < 89 + IMG_WIDTH * IMG_HEIGHT * 6) {
				//copy: blue color frame 0 + green color frame 0 + red color frame 0 + blue color frame 1 + ...
				// color base + color offset + pixel offset (img offset Implicitly through color offset)
				uint32_t copy_length = (byte_left < 89 + IMG_WIDTH * IMG_HEIGHT * 6 - rec_index) ? byte_left : 89 + IMG_WIDTH * IMG_HEIGHT * 6 - rec_index;
				uint8_t *pixel_address = shared_memory->stereo_frame[next_frame_id(shared_memory)].left_blue + rec_index - 89; // IMG_WIDTH * IMG_HEIGHT * index + rec_index - IMG_WIDTH * IMG_HEIGHT * index;
				memcpy(pixel_address, &payload[i], copy_length);
				
				i += copy_length;
				rec_index += copy_length;
			} else {
				xil_printf("WARING: invalid TCP data len error\n\r");
			}
		}
    }
	if(rec_index == 89 + IMG_WIDTH * IMG_HEIGHT * 6){
		shared_memory->stereo_frame[next_frame_id(shared_memory)].status=READY;
		xil_printf("INFO: receive image: %d\n\r", shared_memory->stereo_frame[next_frame_id(shared_memory)].frame_id);
		rec_index = 1;
		connection_status = WAIT;
	}
}

static err_t recv_callback(void *arg, struct tcp_pcb *tpcb,
                               struct pbuf *p, err_t err)
{
	/* do not read the packet if we are not in ESTABLISHED state */
	if (!p) {
		tcp_close(tpcb);
		tcp_recv(tpcb, NULL);
		return ERR_OK;
	}

	/* indicate that the packet has been received */
	tcp_recved(tpcb, p->len);

	if (connection_status == WAIT || connection_status == WAIT_FOR_FRAME) {
		if (p->len > 0){
			uint8_t type = *(uint8_t*)p->payload;
			if(type == 1 || type == 2){
				if(connection_status != WAIT_FOR_FRAME){
					xil_printf("WARNING: receive frame in wrong state: %d", connection_status);
				}
				connection_status = (type == 1) ? RECEIVE_1 : RECEIVE_2;
				pbuf_header(p, -1);   // remove first (read) byte
				recv_img(arg, tpcb, p, err);
			} else if (type == 3) {
				// confirm receipt
				// not used at moment
			} else if (type == 0) {
				xil_printf("INFO: receive notification to close connection (probably due to the completed test run)\n\r");
				connection_status = NOTIFY_TO_CLOSE;
			} else{
				xil_printf("WARNING: receive illegal type: %d\n\r", type);
			}
		}else {
			xil_printf("WARNING: message len < 1\n\r");
		}
	} else if (connection_status == RECEIVE_1 || connection_status == RECEIVE_2) {
		recv_img(arg, tpcb, p, err);
	} else if(connection_status == CLOSE) {
		tcp_client_close(tpcb);
		xil_printf("WARNING: receive message in CLOSE status\n\r");
	} else {
		// illegal send time
		xil_printf("WARNING: receive img at invalid time and status: %d\n\r", connection_status);
	}

	/* free the received pbuf */
	pbuf_free(p);
	return ERR_OK;
}

/** TCP connected callback (active connection), send data now */
static err_t tcp_client_connected(void *arg, struct tcp_pcb *tpcb, err_t err)
{
	if (err != ERR_OK) {
		tcp_client_close(tpcb);
		xil_printf("ERROR: Connection error\n\r");
		return err;
	}
	/* store state */
	c_pcb = tpcb;

	/* set callback values & functions */
	tcp_arg(c_pcb, NULL);
	tcp_sent(c_pcb, tcp_client_sent);
	tcp_err(c_pcb, tcp_client_err);
	tcp_recv(tpcb, recv_callback);

	/* initiate data transfer */
	connection_status = WAIT;
	xil_printf("INFO: Successfully connected to server\r\n");
	return ERR_OK;
}

void start_tcp_client(SharedMemory_t *frame_buffer)
{
	err_t err;
	struct tcp_pcb* pcb;
	ip_addr_t remote_addr;

	shared_memory = frame_buffer;
	connection_status = CLOSE;
#if LWIP_IPV6==1
	remote_addr.type= IPADDR_TYPE_V6;
	err = inet6_aton(TCP_SERVER_IPV6_ADDRESS, &remote_addr);
#else
	err = inet_aton(TCP_SERVER_IP_ADDRESS, &remote_addr);
#endif /* LWIP_IPV6 */

	if (!err) {
		xil_printf("ERROR: Invalid Server IP address: %d\r\n", err);
		return;
	}

	/* Create Client PCB */
	pcb = tcp_new_ip_type(IPADDR_TYPE_ANY);
	if (!pcb) {
		xil_printf("ERROR: Error in PCB creation. out of memory\r\n");
		return;
	}
	tcp_err(pcb, tcp_client_err);
	xil_printf("INFO: try to connect to: %s\r\n", TCP_SERVER_IP_ADDRESS);
	err = tcp_connect(pcb, &remote_addr, TCP_CONN_PORT,
			tcp_client_connected);
	if (err) {
		xil_printf("ERROR: Error on tcp_connect: %d\r\n", err);
		tcp_client_close(pcb);
		return;
	}

	return;
}
