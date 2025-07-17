# ── 컴파일러 및 옵션 ─────────────────────────────────────
CXX       := g++
CXXFLAGS  := -O2 -std=c++17 -static-libgcc -static-libstdc++

# ── 소스, 오브젝트, 디펜던시 ───────────────────────────────
SRC    := main.cpp \
          Momentum.cpp \
          Poisson.cpp \
          TDMASystem.cpp

OBJ    := $(SRC:.cpp=.o)
DEP    := $(SRC:.cpp=.d)
TARGET := main.exe

.PHONY: all clean run mkdir

# ── 기본 빌드 ─────────────────────────────────────────────
all: mkdir $(TARGET)

$(TARGET): $(OBJ)
	@echo "[Link] $@"
	$(CXX) $(CXXFLAGS) -o $@ $^

# ── 패턴 규칙 + 자동 의존성 생성 ───────────────────────────
%.o: %.cpp
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEP)

# ── 실행 ───────────────────────────────────────────────────
run: $(TARGET)
	@echo "[Run] ./$(TARGET)"
	./$(TARGET)

# ── 디렉토리 생성 (예: Result/) ─────────────────────────────
mkdir:
	@mkdir -p Result

# ── 정리 ───────────────────────────────────────────────────
clean:
	@echo "[Clean] Removing objects and dependencies"
	rm -f $(TARGET) $(OBJ) $(DEP)
