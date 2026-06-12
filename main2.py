import euler_trail as ff
import parser as ld

result = open('placement_result.txt', mode = 'w')

graph_info = []
graph_info = ld.load_cells('Nangate_15nm.sp')            ####리스트에 연결 정보 넣기

ld.build_adjacency_matrices(graph_info)                  ######pfet, nfet 각각 인접행렬 구성완료

count_degree = [[[] for _ in range(2)] for _ in range(len(graph_info))]
for i in range(len(graph_info)):
    print('start ',graph_info[i].name)
    pstart = []
    nstart = []
    count_degree[i][0] = ff.vertex_degree(len(graph_info[i].p_vertex), graph_info[i].pmatrix)   #count_degree[i][0] : pfet의 degree 정보
    count_degree[i][1] = ff.vertex_degree(len(graph_info[i].n_vertex), graph_info[i].nmatrix)   #count_degree[i][0] : pfet의 degree 정보

    #print('p_degree, n_degree :', count_degree[i][0], count_degree[i][1])

    for length in range(len(count_degree[i][0])):
        if count_degree[i][0][length] % 2 != 0:
            pstart.append(length)
        if length == len(count_degree[i][0]) - 1 and len(pstart) == 0:
            for num in range(len(count_degree[i][0])):
                pstart.append(num)

    for length in range(len(count_degree[i][1])):
        if count_degree[i][1][length] % 2 != 0:
            nstart.append(length)
        if length == len(count_degree[i][1]) - 1 and len(nstart) == 0:
            for num in range(len(count_degree[i][1])):
                nstart.append(num)

######################시작점 찾기 알고리즘 완성######################

    vertex_num_list = []
    vertex_num_list.append([])
    for cnt_edge in range(len(graph_info[i].connect_list[0])):
        pconnection_info = []
        temp_edge = graph_info[i].edge
        temp_pvertex = graph_info[i].p_vertex

        pconnection_info.append(temp_edge.index(graph_info[i].connect_list[0][cnt_edge][0]))
        pconnection_info.append(temp_pvertex.index(graph_info[i].connect_list[0][cnt_edge][1]))
        pconnection_info.append(temp_pvertex.index(graph_info[i].connect_list[0][cnt_edge][2]))
        vertex_num_list[0].append(pconnection_info)

    vertex_num_list.append([])
    for cnt_edge in range(len(graph_info[i].connect_list[1])):
        nconnection_info = []
        temp_edge = graph_info[i].edge
        temp_nvertex = graph_info[i].n_vertex

        nconnection_info.append(temp_edge.index(graph_info[i].connect_list[1][cnt_edge][0]))
        nconnection_info.append(temp_nvertex.index(graph_info[i].connect_list[1][cnt_edge][1]))
        nconnection_info.append(temp_nvertex.index(graph_info[i].connect_list[1][cnt_edge][2]))
        vertex_num_list[1].append(nconnection_info)

##################edge, vertex 숫자화 완료####################

    odd_pnum = 0
    odd_nnum = 0
    for index in range(len(count_degree[i][0])):
        if count_degree[i][0][index] % 2 != 0:
            odd_pnum += 1
    for index in range(len(count_degree[i][1])):
        if count_degree[i][1][index] % 2 != 0:
            odd_nnum += 1

    total_vertex_num = [[], []]
    same_pconnect_list = []
    ans_location = []

    if (odd_nnum > odd_pnum) or (graph_info[i].name == 'NOR4_X2'):
        copy_pedge = ff.ans_oddnum('pfet', odd_pnum, pstart, graph_info[i], vertex_num_list, total_vertex_num)
    elif odd_pnum >= odd_nnum:
        copy_nedge = ff.ans_oddnum('nfet', odd_nnum, nstart, graph_info[i], vertex_num_list, total_vertex_num)

    if (odd_nnum > odd_pnum) or (graph_info[i].name == 'NOR4_X2'):
        for num, ans in enumerate(graph_info[i].ans_p):
            length = len(vertex_num_list[1])
            include, ans_nvertex, dummy_loc, ans_loc = ff.Is_include(graph_info[i].name, ans, vertex_num_list[1], 0, odd_nnum, 0, len(graph_info[i].edge))
            if include == 1:
                ans_edge = ans

                ff.write_placement(ans_edge, graph_info[i].visited[0][num], ans_nvertex, vertex_num_list[1], graph_info[i], result, dummy_loc, ans_loc)

                if len(vertex_num_list[1]) != length:
                    del vertex_num_list[1][-1]

            else:
                continue

    elif odd_pnum >= odd_nnum:
        for num, ans in enumerate(graph_info[i].ans_n):
            length = len(vertex_num_list[0])
            include, ans_pvertex, dummy_loc, ans_loc = ff.Is_include(graph_info[i].name, ans, vertex_num_list[0], 0, odd_pnum, 0, len(graph_info[i].edge))
            if include == 1:
                ans_edge = ans

                ff.write_placement(ans_edge, ans_pvertex, graph_info[i].visited[1][num], vertex_num_list[0], graph_info[i], result, dummy_loc, ans_loc)

                if len(vertex_num_list[0]) != length:
                    del vertex_num_list[0][-1]

            else:
                continue

    print('end ',graph_info[i].name)

result.close()

result = open('routing_pin_result.txt', mode = 'w')

file = open("placement_result.txt", "r")

lines = file.readlines()
pmos = []
nmos = []
gate = []
for data in lines:
    if ('Name' in data):
        tmp_data = data.strip().split(' ')
        result.write('Name :')
        result.write(" %s" % tmp_data[2])
        result.write("\n")
    elif ('PMOS' in data):
        tmp_data = data.strip().split(' ')[2:]
        pmos = data.strip().split(' ')[2:]
        p_cnt = len(tmp_data)

    elif ('NMOS' in data):
        tmp_data = data.strip().split(' ')[2:]
        nmos = data.strip().split(' ')[2:]
        n_cnt = len(tmp_data)

    elif ('Gate' in data):
        tmp_data = data.strip().split(' ')[2:]
        gate = data.strip().split(' ')[2:]
        gate_cnt = len(tmp_data)
        result.write("%s " %(gate_cnt + p_cnt))
        result.write("10 2 ")
        result.write("\n")

        routing_list = ff.routing_all(pmos, nmos, gate)
        result.write("%s " %(len(routing_list)))
        #print(routing_list)
        result.write("\n")

        ff.write_routing(routing_list, result)

        result.write("\n")
result.close()

cnt = 0
for i in range(len(graph_info)):
    print('')
    print(f"name , degree, odd_pnum, odd_nnum: {graph_info[i].name}, {count_degree[i]}, {odd_pnum}, {odd_nnum}")
    cnt += 1
print('print되는 gate 수 : ', cnt)
