################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../dircrawl.cpp \
../dumpcrawl.cpp \
../fasta.cpp \
../main.cpp \
../parselin.cpp \
../progress.cpp 

OBJS += \
./dircrawl.o \
./dumpcrawl.o \
./fasta.o \
./main.o \
./parselin.o \
./progress.o 

CPP_DEPS += \
./dircrawl.d \
./dumpcrawl.d \
./fasta.d \
./main.d \
./parselin.d \
./progress.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


