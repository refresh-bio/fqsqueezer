all: fqs-0.1

FQS_ROOT_DIR = .
FQS_MAIN_DIR = fqs
FQS_LIBS_DIR = .

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++14 -pthread -mavx -I $(FQS_LIBS_DIR)
CLINK	= -lm -O3 -std=c++14 -pthread -mavx -fabi-version=6 

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

fqs-1.1: $(FQS_MAIN_DIR)/application.o \
	$(FQS_MAIN_DIR)/code_ctx.o \
	$(FQS_MAIN_DIR)/compressor.o \
	$(FQS_MAIN_DIR)/dna.o \
	$(FQS_MAIN_DIR)/ht_kmer.o \
	$(FQS_MAIN_DIR)/id.o \
	$(FQS_MAIN_DIR)/fqsqueezer.o \
	$(FQS_MAIN_DIR)/meta.o \
	$(FQS_MAIN_DIR)/mtf.o \
	$(FQS_MAIN_DIR)/quality.o \
	$(FQS_MAIN_DIR)/utils.o
	$(CC) $(CLINK) -o $(FQS_ROOT_DIR)/$@  \
	$(FQS_MAIN_DIR)/application.o \
	$(FQS_MAIN_DIR)/code_ctx.o \
	$(FQS_MAIN_DIR)/compressor.o \
	$(FQS_MAIN_DIR)/dna.o \
	$(FQS_MAIN_DIR)/ht_kmer.o \
	$(FQS_MAIN_DIR)/id.o \
	$(FQS_MAIN_DIR)/fqsqueezer.o \
	$(FQS_MAIN_DIR)/meta.o \
	$(FQS_MAIN_DIR)/mtf.o \
	$(FQS_MAIN_DIR)/quality.o \
	$(FQS_MAIN_DIR)/utils.o

clean:
	-rm $(FQS_MAIN_DIR)/*.o
	-rm $(FQS_LIBS_DIR)/*.o
	-rm fqs-1.1
	
