#ifndef VECTOR_H
#define VECTOR_H

#include <cstdint>
#include <cstddef> // for size_t
#include <utility> // std::pair

template <typename T, size_t MAX_SIZE>
struct StaticVector {
    T data[MAX_SIZE];
    size_t sz = 0;

    // Anzahl Elemente
    size_t size() const { return sz; }

    // Prüfen ob leer
    bool empty() const { return sz == 0; }

    // Leeren
    void clear() { sz = 0; }

    // Element hinzufügen
    void push_back(const T& value) {
        if(sz < MAX_SIZE) {
            data[sz++] = value;
        }
        // Optional: else Fehlerbehandlung
    }

    // Emplace_back für Konstruktion vor Ort (z.B. std::pair)
    template <typename... Args>
    void emplace_back(Args&&... args) {
        if(sz < MAX_SIZE) {
            data[sz++] = T(std::forward<Args>(args)...);
        }
        // Optional: else Fehlerbehandlung
    }

    // Zugriff auf Elemente
    T& operator[](size_t idx) { return data[idx]; }
    const T& operator[](size_t idx) const { return data[idx]; }

    // Iteratoren
    T* begin() { return data; }
    T* end() { return data + sz; }
};


template <typename T, size_t MAX_SIZE>
struct RingBuffer {
    T ring_buf[MAX_SIZE];
    size_t q_head = 0;
    size_t q_tail = 0;
    size_t q_size = 0;
    size_t index = 0;

    void clear() {
        q_head = q_tail = q_size = index = 0;
    }

    bool is_empty() { return q_size == 0; }
    size_t size() { return q_size; }
    void push(T v) {
        ring_buf[q_tail] = v;
        q_tail++;
        if (q_tail >= MAX_SIZE){
            q_tail = 0;
            if (q_head == 0) {
                q_head = 1; // overwrite oldest
            }
        }
        if (q_tail == q_head) {
            q_head++;
            if (q_head >= MAX_SIZE) {
                q_head = 0;
            }
        }
        q_size++;
    }

    // don't use next after pop call!
    /* T pop() {
        T v = ring_buf[q_head];
        q_head++;
        if (q_head >= MAX_SIZE) q_head = 0;
        q_size--;
        return v;
    }*/
    T next() {
        T v = ring_buf[index];
        index++;
        if (index >= MAX_SIZE) index = 0;
        return v;
    }
    bool has_next() {
        return index != q_tail;
    }
    void reset_index() {
        index = q_head;
    }

};

#endif // VECTOR_H