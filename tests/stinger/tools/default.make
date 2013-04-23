#TARGETS AND FIXED CODE
OBJ=$(addprefix $(OBJDIR)/,$(addsuffix .o,$(basename $(CSRC))))
CFLAGS+= $(addprefix -I,$(INCLUDES))

.PHONY: all clean
all: $(APPNAME)

clean: 
	rm -rf $(APPNAME) $(OBJ)

$(APPNAME): $(OBJ)
ifeq ($(IS_LIB),)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
endif
ifeq ($(IS_LIB),0)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
endif
ifeq ($(IS_LIB),no)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
endif
	
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $^ -o $@ $(LIBS)
