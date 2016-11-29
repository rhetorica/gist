################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../fgs/hmm_lib.cpp \
../fgs/run.cpp \
../fgs/util_lib.cpp 

OBJS += \
./fgs/hmm_lib.o \
./fgs/run.o \
./fgs/util_lib.o 

CPP_DEPS += \
./fgs/hmm_lib.d \
./fgs/run.d \
./fgs/util_lib.d 


# Each subdirectory must supply rules for building sources it contributes
fgs/%.o: ../fgs/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


