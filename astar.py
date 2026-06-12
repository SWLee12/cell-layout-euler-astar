import heapq

def manhattan_distance(a, b):
    x1, y1 = a
    x2, y2 = b
    distance = abs(x2 - x1) + abs(y2 - y1)
    return distance


def reconstruct_path(parent, current):
    path = []
    while current is not None:
        path.append(current)
        current = parent[current]
    return path[::-1]  # 뒤집어서 출발→도착 순으로


def astar(grid, start, goal):
    # grid  : 2D 리스트, 0 = 빈 칸, 1 = 장애물
    # start : 출발 좌표 (x, y)
    # goal  : 목적지 좌표 (x, y)
    
    # open_list = [(f=0, start)] ← 탐색할 셀 목록 (우선순위 큐, f 낮은 것 먼저)
    # closed_set = {} ← 이미 탐색 완료한 셀 목록
    # g_score = {start: 0} ← 출발점→각 셀까지 실제 비용
    # parent = {start: None} ← 경로 역추적용 (어디서 왔는지)
    
    open_list  = []
    closed_set = set()
    g_score    = {start: 0}
    parent     = {start: None}
    
    heapq.heappush(open_list, (manhattan_distance(start, goal), start))

    rows = len(grid)
    cols = len(grid[0])

    while open_list:
        f, current = heapq.heappop(open_list)
        if current == goal:
            return reconstruct_path(parent, current)
        closed_set.add(current)
        
        directions = [(0,1), (0,-1), (1,0), (-1,0)]
        
        for dx, dy in directions:
            neighbor = (current[0] + dx, current[1] + dy)   # 1. 먼저 계산

            if not (0 <= neighbor[0] < cols and 0 <= neighbor[1] < rows):  # 2. 범위 체크
                continue
            if neighbor in closed_set:                        # 3. closed 체크
                continue
            if grid[neighbor[1]][neighbor[0]] == 1:           # 4. 장애물 체크
                continue
            
            new_g = g_score[current] + 1                      # 현재 셀까지의 g + 한 칸 이동 비용(1)
            
            if new_g < g_score.get(neighbor, float('inf')):   # 이 neighbor를 처음 보거나, 더 좋은 경로를 찾았으면
                
                g_score[neighbor] = new_g                     # g_score 업데이트
                
                h = manhattan_distance(neighbor, goal)        # h 계산 (neighbor → goal 거리)
                
                f = new_g + h                                 # f 계산
                
                parent[neighbor] = current                    # parent 기록
                
                # open_list에 추가
                heapq.heappush(open_list, (f, neighbor))
    return None


def parse_routing_file(filename):
    cells = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('Name'):
            cell_name = line.split(':')[1].strip()
            i += 1
            i += 1  # 핀 수 줄 건너뜀
            net_count = int(lines[i].strip())
            i += 1
            
            nets = []
            for _ in range(net_count):
                nums = list(map(int, lines[i].strip().split()))                    # 한 줄에 있는 숫자들을 3개씩 묶어서 (x, y, layer) 로
                pins = [(nums[j], nums[j+1]) for j in range(0, len(nums), 3)]
                nets.append(pins)
                i += 1
            
            cells.append({'name': cell_name, 'nets': nets})
        else:
            i += 1
    
    return cells


def connect_nets(cells):
    max_x, max_y = 0, 0
    for cell in cells:
        for pins in cell['nets']:
            for x, y in pins:
                max_x = max(max_x, x)
                max_y = max(max_y, y)

    for cell in cells:
        grid = [[0] * (max_x + 5) for _ in range(max_y + 5)]
        print(f"\n=== {cell['name']} ===")
        for net_idx, pins in enumerate(cell['nets']):
            print(f"net {net_idx+1}: {pins}")
            net_path = []  # 이 net의 전체 경로 저장
            
            for i in range(len(pins) - 1):
                start = pins[i]
                goal  = pins[i+1]
                path  = astar(grid, start, goal)
                print(f"  {start} → {goal} : {path}")
                
                if path:
                    net_path.extend(path)
            
            # net 전체 연결 끝난 후에 장애물 등록
            for x, y in net_path:
                grid[y][x] = 1

cells = parse_routing_file('routing_pin_result.txt')
connect_nets(cells)    
