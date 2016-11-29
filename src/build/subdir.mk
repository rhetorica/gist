################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Class.cpp \
../DifferenceTable.cpp \
../SparseTable.cpp \
../Table.cpp \
../autocross.cpp \
../build.cpp \
../bwa.cpp \
../classify.cpp \
../main.cpp \
../output.cpp \
../score.cpp 

OBJS += \
./Class.o \
./DifferenceTable.o \
./SparseTable.o \
./Table.o \
./autocross.o \
./build.o \
./bwa.o \
./classify.o \
./main.o \
./output.o \
./score.o 

CPP_DEPS += \
./Class.d \
./DifferenceTable.d \
./SparseTable.d \
./Table.d \
./autocross.d \
./build.d \
./bwa.d \
./classify.d \
./main.d \
./output.d \
./score.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


