#define check_tsk_error(val)                                                         \
    do {                                                                             \
        if (val < 0)                                                                 \
            Rf_error("file %s, line %d: %s", __FILE__, __LINE__, tsk_strerror(val)); \
    } while (0)

#define TSX_ERROR(msg) Rf_error("file %s, line %d: %s", __FILE__, __LINE__, (msg))

#define TSX_WARN(msg) Rf_warning("file %s, line %d: %s", __FILE__, __LINE__, (msg))
    