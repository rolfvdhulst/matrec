#include "matrec/Graphic.h"
#include <assert.h>

//Columns 0..x correspond to elements 0..x
//Rows 0..y correspond to elements -1.. -y-1
#define MARKER_ROW_ELEMENT (INT_MIN)
#define MARKER_COLUMN_ELEMENT (INT_MAX)
typedef int spqr_element;

static bool SPQRelementIsRow(spqr_element element){
    return element < 0;
}
static bool SPQRelementIsColumn(spqr_element element){
    return !SPQRelementIsRow(element);
}
static MATREC_row SPQRelementToRow(spqr_element element){
    assert(SPQRelementIsRow(element));
    return (MATREC_row) (-element - 1);
}
static MATREC_col SPQRelementToColumn(spqr_element element){
    assert(SPQRelementIsColumn(element));
    return (MATREC_col) element;
}

static spqr_element MATRECrowToElement(MATREC_row row){
    assert(MATRECrowIsValid(row));
    return (spqr_element) (-row-1);
}
static spqr_element MATRECcolumnToElement(MATREC_col column){
    assert(MATRECcolIsValid(column));
    return (spqr_element) column;
}

typedef int spqr_node;
#define SPQR_INVALID_NODE (-1)

static bool SPQRnodeIsInvalid(spqr_node node){
    return node < 0;
}
static bool SPQRnodeIsValid(spqr_node node){
    return !SPQRnodeIsInvalid(node);
}

typedef int spqr_member;
#define SPQR_INVALID_MEMBER (-1)

static bool SPQRmemberIsInvalid(spqr_member member){
    return member < 0;
}
static bool SPQRmemberIsValid(spqr_member member){
    return !SPQRmemberIsInvalid(member);
}

typedef int spqr_edge;
#define SPQR_INVALID_EDGE (INT_MAX)

static bool SPQRedgeIsInvalid(spqr_edge edge){
    return edge == SPQR_INVALID_EDGE;
}
static bool SPQRedgeIsValid(spqr_edge edge){
    return !SPQRedgeIsInvalid(edge);
}


typedef enum {
    SPQR_MEMBERTYPE_RIGID = 0, //Also known as triconnected components
    SPQR_MEMBERTYPE_PARALLEL = 1,//Also known as a 'bond'
    SPQR_MEMBERTYPE_SERIES = 2, //Also known as 'polygon' or 'cycle'
    SPQR_MEMBERTYPE_LOOP = 3,
    SPQR_MEMBERTYPE_UNASSIGNED = 4 // To indicate that the member has been merged/is not representative; this is just there to catch errors.
} SPQRMemberType;

typedef struct{
    spqr_edge previous;
    spqr_edge next;
}SPQRGraphicDecompositionEdgeListNode;

typedef struct {
    spqr_node representativeNode;
    spqr_edge firstEdge;//first edge of the neighbouring edges
    int numEdges;
} SPQRGraphicDecompositionNode;

typedef struct {
    spqr_node head;
    spqr_node tail;
    spqr_member member;
    spqr_member childMember;
    SPQRGraphicDecompositionEdgeListNode headEdgeListNode;
    SPQRGraphicDecompositionEdgeListNode tailEdgeListNode;
    SPQRGraphicDecompositionEdgeListNode edgeListNode; //Linked-list node of the array of edges of the member which this edge is in

    spqr_element element;
} SPQRGraphicDecompositionEdge;


typedef struct {
    spqr_member representativeMember;
    SPQRMemberType type;

    spqr_member parentMember;
    spqr_edge markerToParent;
    spqr_edge markerOfParent;

    spqr_edge firstEdge; //First of the members' linked-list edge array
    int num_edges;
} SPQRGraphicDecompositionMember;

struct MATRECGraphicDecompositionImpl {
    int numEdges;
    int memEdges;
    SPQRGraphicDecompositionEdge *edges;
    spqr_edge firstFreeEdge;

    int memMembers;
    int numMembers;
    SPQRGraphicDecompositionMember *members;

    int memNodes;
    int numNodes;
    SPQRGraphicDecompositionNode *nodes;

    int memRows;
    int numRows;
    spqr_edge * rowEdges;

    int memColumns;
    int numColumns;
    spqr_edge * columnEdges;

    MATREC * env;

    int numConnectedComponents;
};

static void swap_ints(int* a, int* b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

static bool nodeIsRepresentative(const MATRECGraphicDecomposition *dec, spqr_node node) {
    assert(dec);
    assert(node < dec->memNodes);
    assert(SPQRnodeIsValid(node));

    return SPQRnodeIsInvalid(dec->nodes[node].representativeNode);
}

static spqr_node findNode(MATRECGraphicDecomposition *dec, spqr_node node) {
    assert(dec);
    assert(SPQRnodeIsValid(node));
    assert(node < dec->memNodes);

    spqr_node current = node;
    spqr_node next;

    //traverse down tree to find the root
    while (SPQRnodeIsValid(next = dec->nodes[current].representativeNode)) {
        current = next;
        assert(current < dec->memNodes);
    }

    spqr_node root = current;
    current = node;

    //update all pointers along path to point to root, flattening the tree
    while (SPQRnodeIsValid(next = dec->nodes[current].representativeNode)) {
        dec->nodes[current].representativeNode = root;
        current = next;
        assert(current < dec->memNodes);
    }
    return root;
}

static spqr_node findNodeNoCompression(const MATRECGraphicDecomposition *dec, spqr_node node) {
    assert(dec);
    assert(SPQRnodeIsValid(node));
    assert(node < dec->memNodes);

    spqr_node current = node;
    spqr_node next;

    //traverse down tree to find the root
    while (SPQRnodeIsValid(next = dec->nodes[current].representativeNode)) {
        current = next;
        assert(current < dec->memNodes);
    }
    spqr_node root = current;
    return root;
}

static spqr_node findEdgeTail(MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_node representative = findNode(dec, dec->edges[edge].tail);
    dec->edges[edge].tail = representative; //update the edge information

    return representative;
}
static spqr_node findEdgeHead(MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_node representative = findNode(dec, dec->edges[edge].head);
    dec->edges[edge].head = representative;//update the edge information

    return representative;
}
static spqr_node findEdgeHeadNoCompression(const MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_node representative = findNodeNoCompression(dec, dec->edges[edge].head);
    return representative;
}
static spqr_node findEdgeTailNoCompression(const MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_node representative = findNodeNoCompression(dec, dec->edges[edge].tail);
    return representative;
}

static spqr_edge getFirstNodeEdge(const MATRECGraphicDecomposition * dec, spqr_node node){
    assert(dec);
    assert(SPQRnodeIsValid(node));
    assert(node < dec->memNodes);
    return dec->nodes[node].firstEdge;
}
static spqr_edge getNextNodeEdgeNoCompression(const MATRECGraphicDecomposition * dec, spqr_edge edge, spqr_node node){
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);
    assert(nodeIsRepresentative(dec,node));

    if(findEdgeHeadNoCompression(dec,edge) == node){
        edge = dec->edges[edge].headEdgeListNode.next;
    }else{
        assert(findEdgeTailNoCompression(dec,edge) == node);
        edge = dec->edges[edge].tailEdgeListNode.next;
    }
    return edge;
}
static spqr_edge getNextNodeEdge(MATRECGraphicDecomposition * dec, spqr_edge edge, spqr_node node){
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);
    assert(nodeIsRepresentative(dec,node));

    if(findEdgeHead(dec,edge) == node){
        edge = dec->edges[edge].headEdgeListNode.next;
    }else{
        assert(findEdgeTailNoCompression(dec,edge) == node);
        dec->edges[edge].tail = node; //This assignment is not necessary but speeds up future queries.
        edge = dec->edges[edge].tailEdgeListNode.next;
    }
    return edge;
}
static spqr_edge getPreviousNodeEdge(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_node node){
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);
    assert(nodeIsRepresentative(dec,node));

    if(findEdgeHead(dec,edge) == node){
        edge = dec->edges[edge].headEdgeListNode.previous;
    }else{
        assert(findEdgeTailNoCompression(dec,edge) == node);
        dec->edges[edge].tail = node; //This assignment is not necessary but speeds up future queries.
        edge = dec->edges[edge].tailEdgeListNode.previous;
    }
    return edge;
}

static void mergeNodeEdgeList(MATRECGraphicDecomposition *dec, spqr_node toMergeInto, spqr_node toRemove){

    spqr_edge firstIntoEdge = getFirstNodeEdge(dec, toMergeInto);
    spqr_edge firstFromEdge = getFirstNodeEdge(dec, toRemove);
    if(SPQRedgeIsInvalid(firstIntoEdge)){
        //new node has no edges
        dec->nodes[toMergeInto].numEdges += dec->nodes[toRemove].numEdges;
        dec->nodes[toRemove].numEdges = 0;

        dec->nodes[toMergeInto].firstEdge = dec->nodes[toRemove].firstEdge;
        dec->nodes[toRemove].firstEdge = SPQR_INVALID_EDGE;

        return;
    }else if (SPQRedgeIsInvalid(firstFromEdge)){
        //Old node has no edges; we can just return
        return;
    }

    spqr_edge lastIntoEdge = getPreviousNodeEdge(dec, firstIntoEdge, toMergeInto);
    assert(SPQRedgeIsValid(lastIntoEdge));
    spqr_edge lastFromEdge = getPreviousNodeEdge(dec, firstFromEdge, toRemove);
    assert(SPQRedgeIsValid(lastFromEdge));


    SPQRGraphicDecompositionEdgeListNode * firstIntoNode = findEdgeHead(dec, firstIntoEdge) == toMergeInto ?
                                                           &dec->edges[firstIntoEdge].headEdgeListNode :
                                                           &dec->edges[firstIntoEdge].tailEdgeListNode;
    SPQRGraphicDecompositionEdgeListNode * lastIntoNode = findEdgeHead(dec, lastIntoEdge) == toMergeInto ?
                                                          &dec->edges[lastIntoEdge].headEdgeListNode :
                                                          &dec->edges[lastIntoEdge].tailEdgeListNode;

    SPQRGraphicDecompositionEdgeListNode * firstFromNode = findEdgeHead(dec, firstFromEdge) == toRemove ?
                                                           &dec->edges[firstFromEdge].headEdgeListNode :
                                                           &dec->edges[firstFromEdge].tailEdgeListNode;
    SPQRGraphicDecompositionEdgeListNode * lastFromNode = findEdgeHead(dec, lastFromEdge) == toRemove ?
                                                          &dec->edges[lastFromEdge].headEdgeListNode :
                                                          &dec->edges[lastFromEdge].tailEdgeListNode;

    firstIntoNode->previous = lastFromEdge;
    lastIntoNode->next = firstFromEdge;
    firstFromNode->previous = lastIntoEdge;
    lastFromNode->next = firstIntoEdge;

    dec->nodes[toMergeInto].numEdges += dec->nodes[toRemove].numEdges;
    dec->nodes[toRemove].numEdges = 0;
    dec->nodes[toRemove].firstEdge = SPQR_INVALID_EDGE;
}

static spqr_node mergeNodes(MATRECGraphicDecomposition *dec, spqr_node first, spqr_node second) {
    assert(dec);
    assert(nodeIsRepresentative(dec, first));
    assert(nodeIsRepresentative(dec, second));
    assert(first != second); //We cannot merge a node into itself
    assert(first < dec->memNodes);
    assert(second < dec->memNodes);

    //The rank is stored as a negative number: we decrement it making the negative number larger.
    // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
    spqr_node firstRank = dec->nodes[first].representativeNode;
    spqr_node secondRank = dec->nodes[second].representativeNode;
    if (firstRank > secondRank) {
        swap_ints(&first, &second);
    }
    //first becomes representative; we merge all of the edges of second into first
    mergeNodeEdgeList(dec,first,second);
    dec->nodes[second].representativeNode = first;
    if (firstRank == secondRank) {
        --dec->nodes[first].representativeNode;
    }
    return first;
}

static bool memberIsRepresentative(const MATRECGraphicDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(member < dec->memMembers);
    assert(SPQRmemberIsValid(member));

    return SPQRmemberIsInvalid(dec->members[member].representativeMember);
}

static spqr_member findMember(MATRECGraphicDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(member < dec->memMembers);
    assert(SPQRmemberIsValid(member));

    spqr_member current = member;
    spqr_member next;

    //traverse down tree to find the root
    while (SPQRmemberIsValid(next = dec->members[current].representativeMember)) {
        current = next;
        assert(current < dec->memMembers);
    }

    spqr_member root = current;
    current = member;

    //update all pointers along path to point to root, flattening the tree
    while (SPQRmemberIsValid(next = dec->members[current].representativeMember)) {
        dec->members[current].representativeMember = root;
        current = next;
        assert(current < dec->memMembers);
    }
    return root;
}

static spqr_member findMemberNoCompression(const MATRECGraphicDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(member < dec->memMembers);
    assert(SPQRmemberIsValid(member));

    spqr_member current = member;
    spqr_member next;

    //traverse down tree to find the root
    while (SPQRmemberIsValid(next = dec->members[current].representativeMember)) {
        current = next;
        assert(current < dec->memMembers);
    }

    spqr_member root = current;
    return root;
}

static spqr_member mergeMembers(MATRECGraphicDecomposition *dec, spqr_member first, spqr_member second) {
    assert(dec);
    assert(memberIsRepresentative(dec, first));
    assert(memberIsRepresentative(dec, second));
    assert(first != second); //We cannot merge a member into itself
    assert(first < dec->memMembers);
    assert(second < dec->memMembers);

    //The rank is stored as a negative number: we decrement it making the negative number larger.
    // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
    spqr_member firstRank = dec->members[first].representativeMember;
    spqr_member secondRank = dec->members[second].representativeMember;
    if (firstRank > secondRank) {
        swap_ints(&first, &second);
    }
    dec->members[second].representativeMember = first;
    if (firstRank == secondRank) {
        --dec->members[first].representativeMember;
    }
    return first;
}





static spqr_member findEdgeMember(MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_member representative = findMember(dec, dec->edges[edge].member);
    dec->edges[edge].member = representative;
    return representative;
}

static spqr_member findEdgeMemberNoCompression(const MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_member representative = findMemberNoCompression(dec, dec->edges[edge].member);
    return representative;
}

static spqr_member findMemberParent(MATRECGraphicDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(member < dec->memMembers);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));


    if(SPQRmemberIsInvalid(dec->members[member].parentMember)){
        return dec->members[member].parentMember;
    }
    spqr_member parent_representative = findMember(dec, dec->members[member].parentMember);
    dec->members[member].parentMember = parent_representative;

    return parent_representative;
}

static spqr_member findMemberParentNoCompression(const MATRECGraphicDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(member < dec->memMembers);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));

    if(SPQRmemberIsInvalid(dec->members[member].parentMember)){
        return dec->members[member].parentMember;
    }
    spqr_member parent_representative = findMemberNoCompression(dec, dec->members[member].parentMember);
    return parent_representative;
}

static spqr_member findEdgeChildMember(MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_member representative = findMember(dec, dec->edges[edge].childMember);
    dec->edges[edge].childMember = representative;
    return representative;
}

static spqr_member findEdgeChildMemberNoCompression(const MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    spqr_member representative = findMemberNoCompression(dec, dec->edges[edge].childMember);
    return representative;
}

//TODO: fix usages, is misleading. Only accounts for CHILD markers, not parent markers!
static bool edgeIsMarker(const MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    return SPQRmemberIsValid(dec->edges[edge].childMember);
}

static bool edgeIsTree(const MATRECGraphicDecomposition *dec, spqr_edge edge) {
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    return SPQRelementIsRow(dec->edges[edge].element);
}

static spqr_element edgeGetElement(const MATRECGraphicDecomposition * dec, spqr_edge edge){
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);

    return dec->edges[edge].element;
}
bool MATRECGraphicDecompositionContainsRow(const MATRECGraphicDecomposition * dec, MATREC_row row){
    assert(MATRECrowIsValid(row) && (int) row < dec->memRows);
    assert(dec);
    return SPQRedgeIsValid(dec->rowEdges[row]);
}
bool MATRECGraphicDecompositionContainsColumn(const MATRECGraphicDecomposition *dec, MATREC_col col){
    assert(MATRECcolIsValid(col) && (int) col < dec->memColumns);
    assert(dec);
    return SPQRedgeIsValid(dec->columnEdges[col]);
}
static void setDecompositionColumnEdge(MATRECGraphicDecomposition *dec, MATREC_col col, spqr_edge edge){
    assert(MATRECcolIsValid(col) && (int)col < dec->memColumns);
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    dec->columnEdges[col] = edge;
}
static void setDecompositionRowEdge(MATRECGraphicDecomposition *dec, MATREC_row row, spqr_edge edge){
    assert(MATRECrowIsValid(row) && (int) row < dec->memRows);
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    dec->rowEdges[row] = edge;
}
static spqr_edge getDecompositionColumnEdge(const MATRECGraphicDecomposition *dec, MATREC_col col){
    assert(MATRECcolIsValid(col) && (int) col < dec->memColumns);
    assert(dec);
    return dec->columnEdges[col];
}
static spqr_edge getDecompositionRowEdge(const MATRECGraphicDecomposition *dec, MATREC_row row){
    assert(MATRECrowIsValid(row) && (int) row < dec->memRows);
    assert(dec);
    return dec->rowEdges[row];
}

MATREC_ERROR MATRECGraphicDecompositionCreate(MATREC * env, MATRECGraphicDecomposition **pDecomposition, int numRows, int numColumns){
    assert(env);
    assert(pDecomposition);
    assert(!*pDecomposition);

    MATREC_CALL(MATRECallocBlock(env, pDecomposition));
    MATRECGraphicDecomposition *dec = *pDecomposition;
    dec->env = env;

    //Initialize edge array data
    int initialMemEdges = 8;
    {
        assert(initialMemEdges > 0);
        dec->memEdges = initialMemEdges;
        dec->numEdges = 0;
        MATREC_CALL(MATRECallocBlockArray(env, &dec->edges, (size_t) dec->memEdges));
        for (spqr_edge i = 0; i < dec->memEdges; ++i) {
            dec->edges[i].edgeListNode.next = i + 1;
            dec->edges[i].member = SPQR_INVALID_MEMBER;
        }
        dec->edges[dec->memEdges - 1].edgeListNode.next = SPQR_INVALID_EDGE;
        dec->firstFreeEdge = 0;
    }

    //Initialize member array data
    int initialMemMembers = 8;
    {
        assert(initialMemMembers > 0);
        dec->memMembers = initialMemMembers;
        dec->numMembers = 0;
        MATREC_CALL(MATRECallocBlockArray(env, &dec->members, (size_t) dec->memMembers));
    }

    //Initialize node array data
    int initialMemNodes = 8;
    {
        assert(initialMemNodes > 0);
        dec->memNodes = initialMemNodes;
        dec->numNodes = 0;
        MATREC_CALL(MATRECallocBlockArray(env, &dec->nodes, (size_t) dec->memNodes));
    }

    //Initialize mappings for rows
    {
        dec->memRows = numRows;
        MATREC_CALL(MATRECallocBlockArray(env, &dec->rowEdges, (size_t) dec->memRows));
        for (int i = 0; i < dec->memRows; ++i) {
            dec->rowEdges[i] = SPQR_INVALID_EDGE;
        }
    }
    //Initialize mappings for columns
    {
        dec->memColumns = numColumns;
        dec->numColumns = 0;
        MATREC_CALL(MATRECallocBlockArray(env, &dec->columnEdges, (size_t) dec->memColumns));
        for (int i = 0; i < dec->memColumns; ++i) {
            dec->columnEdges[i] = SPQR_INVALID_EDGE;
        }
    }

    dec->numConnectedComponents = 0;
    return MATREC_OKAY;
}

void MATRECGraphicDecompositionFree(MATRECGraphicDecomposition **pDec){
    assert(pDec);
    assert(*pDec);

    MATRECGraphicDecomposition *dec = *pDec;
    MATRECfreeBlockArray(dec->env, &dec->columnEdges);
    MATRECfreeBlockArray(dec->env, &dec->rowEdges);
    MATRECfreeBlockArray(dec->env, &dec->nodes);
    MATRECfreeBlockArray(dec->env, &dec->members);
    MATRECfreeBlockArray(dec->env, &dec->edges);

    MATRECfreeBlock(dec->env, pDec);

}
static spqr_edge getFirstMemberEdge(const MATRECGraphicDecomposition * dec, spqr_member member){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    return dec->members[member].firstEdge;
}
static spqr_edge getNextMemberEdge(const MATRECGraphicDecomposition * dec, spqr_edge edge){
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);
    edge = dec->edges[edge].edgeListNode.next;
    return edge;
}
static spqr_edge getPreviousMemberEdge(const MATRECGraphicDecomposition *dec, spqr_edge edge){
    assert(dec);
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);
    edge = dec->edges[edge].edgeListNode.previous;
    return edge;
}

static void addEdgeToMemberEdgeList(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_member member){
    spqr_edge firstMemberEdge = getFirstMemberEdge(dec, member);

    if(SPQRedgeIsValid(firstMemberEdge)){
        spqr_edge lastMemberEdge = getPreviousMemberEdge(dec, firstMemberEdge);
        dec->edges[edge].edgeListNode.next = firstMemberEdge;
        dec->edges[edge].edgeListNode.previous = lastMemberEdge;
        dec->edges[firstMemberEdge].edgeListNode.previous = edge;
        dec->edges[lastMemberEdge].edgeListNode.next = edge;
    }else{
        assert(dec->members[member].num_edges == 0);
        dec->edges[edge].edgeListNode.next = edge;
        dec->edges[edge].edgeListNode.previous = edge;
    }
    dec->members[member].firstEdge = edge;//TODO: update this in case of row/column edges to make memory ordering nicer?
    ++(dec->members[member].num_edges);
}
static MATREC_ERROR createEdge(MATRECGraphicDecomposition *dec, spqr_member member, spqr_edge *pEdge) {
    assert(dec);
    assert(pEdge);
    assert(SPQRmemberIsInvalid(member) || memberIsRepresentative(dec, member));

    spqr_edge index = dec->firstFreeEdge;
    if (SPQRedgeIsValid(index)) {
        dec->firstFreeEdge = dec->edges[index].edgeListNode.next;
    } else {
        //Enlarge array, no free nodes in edge list
        int newSize = 2 * dec->memEdges;
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &dec->edges, (size_t) newSize));
        for (int i = dec->memEdges + 1; i < newSize; ++i) {
            dec->edges[i].edgeListNode.next = i + 1;
            dec->edges[i].member = SPQR_INVALID_MEMBER;
        }
        dec->edges[newSize - 1].edgeListNode.next = SPQR_INVALID_EDGE;
        dec->firstFreeEdge = dec->memEdges + 1;
        index = dec->memEdges;
        dec->memEdges = newSize;
    }
    //TODO: Is defaulting these here necessary?
    dec->edges[index].tail = SPQR_INVALID_NODE;
    dec->edges[index].head = SPQR_INVALID_NODE;
    dec->edges[index].member = member;
    dec->edges[index].childMember = SPQR_INVALID_MEMBER;

    dec->edges[index].headEdgeListNode.next = SPQR_INVALID_EDGE;
    dec->edges[index].headEdgeListNode.previous = SPQR_INVALID_EDGE;
    dec->edges[index].tailEdgeListNode.next = SPQR_INVALID_EDGE;
    dec->edges[index].tailEdgeListNode.previous = SPQR_INVALID_EDGE;

    dec->numEdges++;

    *pEdge = index;

    return MATREC_OKAY;
}
static MATREC_ERROR createRowEdge(MATRECGraphicDecomposition *dec, spqr_member member, spqr_edge *pEdge, MATREC_row row){
    MATREC_CALL(createEdge(dec,member,pEdge));
    setDecompositionRowEdge(dec,row,*pEdge);
    addEdgeToMemberEdgeList(dec,*pEdge,member);
    dec->edges[*pEdge].element = MATRECrowToElement(row);

    return MATREC_OKAY;
}
static MATREC_ERROR createColumnEdge(MATRECGraphicDecomposition *dec, spqr_member member, spqr_edge *pEdge, MATREC_col column){
    MATREC_CALL(createEdge(dec,member,pEdge));
    setDecompositionColumnEdge(dec,column,*pEdge);
    addEdgeToMemberEdgeList(dec,*pEdge,member);
    dec->edges[*pEdge].element = MATRECcolumnToElement(column);

    return MATREC_OKAY;
}
static MATREC_ERROR createMember(MATRECGraphicDecomposition *dec, SPQRMemberType type, spqr_member * pMember){
    assert(dec);
    assert(pMember);

    if(dec->numMembers == dec->memMembers){
        dec->memMembers *= 2;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&dec->members,(size_t) dec->memMembers));
    }
    SPQRGraphicDecompositionMember *data = &dec->members[dec->numMembers];
    data->markerOfParent = SPQR_INVALID_EDGE;
    data->markerToParent = SPQR_INVALID_EDGE;
    data->firstEdge = SPQR_INVALID_EDGE;
    data->representativeMember = SPQR_INVALID_MEMBER;
    data->num_edges = 0;
    data->parentMember = SPQR_INVALID_MEMBER;
    data->type = type;

    *pMember = dec->numMembers;

    dec->numMembers++;
    return MATREC_OKAY;
}

static MATREC_ERROR createNode(MATRECGraphicDecomposition *dec, spqr_node * pNode){

    if(dec->numNodes == dec->memNodes){
        dec->memNodes*=2;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&dec->nodes,(size_t) dec->memNodes));
    }
    *pNode = dec->numNodes;
    dec->nodes[dec->numNodes].representativeNode = SPQR_INVALID_NODE;
    dec->nodes[dec->numNodes].firstEdge = SPQR_INVALID_EDGE;
    dec->nodes[dec->numNodes].numEdges = 0;
    dec->numNodes++;

    return MATREC_OKAY;
}
static void removeEdgeFromNodeEdgeList(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_node node, bool nodeIsHead){
    SPQRGraphicDecompositionEdgeListNode * edgeListNode = nodeIsHead ? &dec->edges[edge].headEdgeListNode : &dec->edges[edge].tailEdgeListNode;

    if(dec->nodes[node].numEdges == 1){
        dec->nodes[node].firstEdge = SPQR_INVALID_EDGE;
    }else{
        spqr_edge next_edge = edgeListNode->next;
        spqr_edge prev_edge = edgeListNode->previous;
        SPQRGraphicDecompositionEdgeListNode * nextListNode = findEdgeHead(dec, next_edge) == node ? &dec->edges[next_edge].headEdgeListNode : &dec->edges[next_edge].tailEdgeListNode;//TODO: finds necessary?
        SPQRGraphicDecompositionEdgeListNode * prevListNode = findEdgeHead(dec, prev_edge) == node ? &dec->edges[prev_edge].headEdgeListNode : &dec->edges[prev_edge].tailEdgeListNode;//TODO: finds necessary?

        nextListNode->previous = prev_edge;
        prevListNode->next = next_edge;

        if(dec->nodes[node].firstEdge == edge){
            dec->nodes[node].firstEdge = next_edge; //TODO: fix this if we want fixed ordering for tree/nontree edges in memory
        }
    }
    //TODO: empty edgeListNode's data? Might not be all that relevant
    --(dec->nodes[node].numEdges);
}
static void addEdgeToNodeEdgeList(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_node node, bool nodeIsHead){
    assert(nodeIsRepresentative(dec,node));

    spqr_edge firstNodeEdge = getFirstNodeEdge(dec, node);

    SPQRGraphicDecompositionEdgeListNode * edgeListNode = nodeIsHead ? &dec->edges[edge].headEdgeListNode : &dec->edges[edge].tailEdgeListNode;
    if(SPQRedgeIsValid(firstNodeEdge)){
        bool nextIsHead = findEdgeHead(dec,firstNodeEdge) == node;
        SPQRGraphicDecompositionEdgeListNode *nextListNode = nextIsHead ? &dec->edges[firstNodeEdge].headEdgeListNode : &dec->edges[firstNodeEdge].tailEdgeListNode;
        spqr_edge lastNodeEdge = nextListNode->previous;

        edgeListNode->next = firstNodeEdge;
        edgeListNode->previous = lastNodeEdge;


        bool previousIsHead = findEdgeHead(dec,lastNodeEdge) == node;
        SPQRGraphicDecompositionEdgeListNode *previousListNode = previousIsHead ? &dec->edges[lastNodeEdge].headEdgeListNode : &dec->edges[lastNodeEdge].tailEdgeListNode;
        previousListNode->next = edge;
        nextListNode->previous = edge;

    }else{
        edgeListNode->next = edge;
        edgeListNode->previous = edge;
    }
    dec->nodes[node].firstEdge = edge; //TODO: update this in case of row/column edges to make memory ordering nicer?er?
    ++dec->nodes[node].numEdges;
    if(nodeIsHead){
        dec->edges[edge].head = node;
    }else{
        dec->edges[edge].tail = node;
    }
}
static void setEdgeHeadAndTail(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_node head, spqr_node tail){
    addEdgeToNodeEdgeList(dec,edge,head,true);
    addEdgeToNodeEdgeList(dec,edge,tail,false);
}
static void clearEdgeHeadAndTail(MATRECGraphicDecomposition *dec, spqr_edge edge){
    removeEdgeFromNodeEdgeList(dec,edge,findEdgeHead(dec,edge),true);
    removeEdgeFromNodeEdgeList(dec,edge,findEdgeTail(dec,edge),false);
    dec->edges[edge].head = SPQR_INVALID_NODE;
    dec->edges[edge].tail = SPQR_INVALID_NODE;
}
static void changeEdgeHead(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_node oldHead, spqr_node newHead){
    assert(nodeIsRepresentative(dec,oldHead));
    assert(nodeIsRepresentative(dec,newHead));
    removeEdgeFromNodeEdgeList(dec,edge,oldHead,true);
    addEdgeToNodeEdgeList(dec,edge,newHead,true);
}
static void changeEdgeTail(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_node oldTail, spqr_node newTail){
    assert(nodeIsRepresentative(dec,oldTail));
    assert(nodeIsRepresentative(dec,newTail));
    removeEdgeFromNodeEdgeList(dec,edge,oldTail,false);
    addEdgeToNodeEdgeList(dec,edge,newTail,false);
}
static void flipEdge(MATRECGraphicDecomposition *dec, spqr_edge edge){
    swap_ints(&dec->edges[edge].head,&dec->edges[edge].tail);

    SPQRGraphicDecompositionEdgeListNode temp = dec->edges[edge].headEdgeListNode;
    dec->edges[edge].headEdgeListNode = dec->edges[edge].tailEdgeListNode;
    dec->edges[edge].tailEdgeListNode = temp;

}
static int nodeDegree(MATRECGraphicDecomposition *dec, spqr_node node){
    assert(dec);
    assert(SPQRnodeIsValid(node));
    assert(node < dec->memNodes);
    return dec->nodes[node].numEdges;
}
static SPQRMemberType getMemberType(const MATRECGraphicDecomposition *dec, spqr_member member){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));
    return dec->members[member].type;
}
static void updateMemberType(const MATRECGraphicDecomposition *dec, spqr_member member, SPQRMemberType type){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));

    dec->members[member].type = type;
}
static spqr_edge markerToParent(const MATRECGraphicDecomposition *dec, spqr_member member){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));
    return dec->members[member].markerToParent;
}
static char typeToChar(SPQRMemberType type){
    switch (type) {
        case SPQR_MEMBERTYPE_RIGID:
            return 'R';
        case SPQR_MEMBERTYPE_PARALLEL:
            return 'P';
        case SPQR_MEMBERTYPE_SERIES:
            return 'S';
        case SPQR_MEMBERTYPE_LOOP:
            return 'L';
        default:
            return '?';
    }
}

static void updateMemberParentInformation(MATRECGraphicDecomposition *dec, const spqr_member newMember, const spqr_member toRemove){
    assert(memberIsRepresentative(dec,newMember));
    assert(findMemberNoCompression(dec,toRemove) == newMember);

    dec->members[newMember].markerOfParent = dec->members[toRemove].markerOfParent;
    dec->members[newMember].markerToParent = dec->members[toRemove].markerToParent;
    dec->members[newMember].parentMember = dec->members[toRemove].parentMember;

    dec->members[toRemove].markerOfParent = SPQR_INVALID_EDGE;
    dec->members[toRemove].markerToParent = SPQR_INVALID_EDGE;
    dec->members[toRemove].parentMember = SPQR_INVALID_MEMBER;
}
static spqr_edge markerOfParent(const MATRECGraphicDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));
    return dec->members[member].markerOfParent;
}



static int getNumMemberEdges(const MATRECGraphicDecomposition * dec, spqr_member member){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));
    return dec->members[member].num_edges;
}

static int getNumNodes(const MATRECGraphicDecomposition *dec){
    assert(dec);
    return dec->numNodes;
}
static int getNumMembers(const MATRECGraphicDecomposition *dec){
    assert(dec);
    return dec->numMembers;
}
static MATREC_ERROR createStandaloneParallel(MATRECGraphicDecomposition *dec, MATREC_col * columns, int num_columns, MATREC_row row, spqr_member * pMember){
    spqr_member member;
    SPQRMemberType type = num_columns < 2 ? SPQR_MEMBERTYPE_LOOP : SPQR_MEMBERTYPE_PARALLEL;
    MATREC_CALL(createMember(dec, type, &member));

    spqr_edge row_edge;
    MATREC_CALL(createRowEdge(dec,member,&row_edge,row));

    spqr_edge col_edge;
    for (int i = 0; i < num_columns; ++i) {
        MATREC_CALL(createColumnEdge(dec,member,&col_edge,columns[i]));
    }
    *pMember = member;

    ++dec->numConnectedComponents;
    return MATREC_OKAY;
}

//TODO: fix tracking connectivity more cleanly, should not be left up to the algorithms ideally
static MATREC_ERROR createConnectedParallel(MATRECGraphicDecomposition *dec, MATREC_col * columns, int num_columns, MATREC_row row, spqr_member * pMember){
    spqr_member member;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &member));

    spqr_edge row_edge;
    MATREC_CALL(createRowEdge(dec,member,&row_edge,row));

    spqr_edge col_edge;
    for (int i = 0; i < num_columns; ++i) {
        MATREC_CALL(createColumnEdge(dec,member,&col_edge,columns[i]));
    }
    *pMember = member;

    return MATREC_OKAY;
}

static MATREC_ERROR createStandaloneSeries(MATRECGraphicDecomposition *dec, MATREC_row * rows, int numRows, MATREC_col col, spqr_member * pMember){
    spqr_member member;
    SPQRMemberType type = numRows < 2 ? SPQR_MEMBERTYPE_LOOP : SPQR_MEMBERTYPE_SERIES;
    MATREC_CALL(createMember(dec, type, &member));

    spqr_edge colEdge;
    MATREC_CALL(createColumnEdge(dec,member,&colEdge,col));

    spqr_edge rowEdge;
    for (int i = 0; i < numRows; ++i) {
        MATREC_CALL(createRowEdge(dec,member,&rowEdge,rows[i]));
    }
    *pMember = member;
    ++dec->numConnectedComponents;
    return MATREC_OKAY;
}
static MATREC_ERROR createConnectedSeries(MATRECGraphicDecomposition *dec, MATREC_row * rows, int numRows, MATREC_col col, spqr_member * pMember){
    spqr_member member;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &member));

    spqr_edge colEdge;
    MATREC_CALL(createColumnEdge(dec,member,&colEdge,col));

    spqr_edge rowEdge;
    for (int i = 0; i < numRows; ++i) {
        MATREC_CALL(createRowEdge(dec,member,&rowEdge,rows[i]));
    }
    *pMember = member;
    return MATREC_OKAY;
}

static void removeEdgeFromMemberEdgeList(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_member member){
    assert(findEdgeMemberNoCompression(dec,edge) == member);
    assert(memberIsRepresentative(dec,member));

    if(dec->members[member].num_edges == 1){
        dec->members[member].firstEdge = SPQR_INVALID_EDGE;

        //TODO: also set edgeListNode to invalid, maybe? Not necessary probably
    }else{
        spqr_edge nextEdge = dec->edges[edge].edgeListNode.next;
        spqr_edge prevEdge = dec->edges[edge].edgeListNode.previous;

        dec->edges[nextEdge].edgeListNode.previous = prevEdge;
        dec->edges[prevEdge].edgeListNode.next = nextEdge;

        if(dec->members[member].firstEdge == edge){
            dec->members[member].firstEdge = nextEdge; //TODO: fix this if we want fixed ordering for tree/nontree edges in memory
        }
    }


    --(dec->members[member].num_edges);
}


static void process_edge(MATREC_row * fundamental_cycle_edges, int * num_cycle_edges, spqr_edge * callStack, int * callStackSize, spqr_edge edge, const MATRECGraphicDecomposition * dec){
    assert(edgeIsTree(dec,edge));
    if(!edgeIsMarker(dec,edge)){
        spqr_member current_member = findEdgeMemberNoCompression(dec, edge);
        if(markerToParent(dec,current_member) == edge){
            spqr_edge other_edge = markerOfParent(dec, current_member);
            assert(!edgeIsTree(dec,other_edge));
            callStack[*callStackSize] = other_edge;
            ++(*callStackSize);
        }else{
            spqr_element element = edgeGetElement(dec,edge);
            assert(SPQRelementIsRow(element));
            fundamental_cycle_edges[*num_cycle_edges] = SPQRelementToRow(element);
            ++(*num_cycle_edges);
        }
    }else{
        spqr_member child_member = findEdgeChildMemberNoCompression(dec, edge);
        spqr_edge other_edge = markerToParent(dec, child_member);
        assert(!edgeIsTree(dec,other_edge));
        callStack[*callStackSize] = other_edge;
        ++(*callStackSize);
    }
}

static int decompositionGetFundamentalCycleRows(const MATRECGraphicDecomposition *dec, MATREC_col column, MATREC_row * output){
    spqr_edge edge = getDecompositionColumnEdge(dec, column);
    if(SPQRedgeIsInvalid(edge)){
        return 0;
    }
    int num_rows = 0;

    spqr_edge * callStack;
    //TODO: probably an overkill amount of memory allocated here... How can we allocate just enough?
    MATREC_ERROR result = MATRECallocBlockArray(dec->env,&callStack,(size_t) dec->memRows);
    if(result != MATREC_OKAY){
        return -1;
    }
    int callStackSize = 1;
    callStack[0] = edge;

    bool * nodeVisited;
    result = MATRECallocBlockArray(dec->env,&nodeVisited,(size_t) dec->numNodes);
    if(result != MATREC_OKAY){
        return -1;
    }
    for (int i = 0; i < dec->numNodes; ++i) {
        nodeVisited[i] = false;
    }

    typedef struct {
        spqr_node node;
        spqr_edge nodeEdge;
    } DFSCallData;
    DFSCallData * pathSearchCallStack;
    result = MATRECallocBlockArray(dec->env,&pathSearchCallStack,(size_t) dec->numNodes);
    if(result != MATREC_OKAY){
        return -1;
    }
    int pathSearchCallStackSize = 0;

    while(callStackSize > 0){
        spqr_edge column_edge = callStack[callStackSize - 1];
        --callStackSize;
        spqr_member column_edge_member = findEdgeMemberNoCompression(dec, column_edge);
        switch(getMemberType(dec,column_edge_member)){
            case SPQR_MEMBERTYPE_RIGID:
            {

                spqr_node source = findEdgeHeadNoCompression(dec, column_edge);
                spqr_node target = findEdgeTailNoCompression(dec, column_edge);

                assert(pathSearchCallStackSize == 0);
                pathSearchCallStack[0].node = source;
                pathSearchCallStack[0].nodeEdge = getFirstNodeEdge(dec,source);
                pathSearchCallStackSize++;
                while(pathSearchCallStackSize > 0){
                    DFSCallData * dfsData  = &pathSearchCallStack[pathSearchCallStackSize-1];
                    nodeVisited[dfsData->node] = true;
                    //cannot be a tree edge which is its parent
                    if(edgeIsTree(dec,dfsData->nodeEdge) &&
                       (pathSearchCallStackSize <= 1 || dfsData->nodeEdge != pathSearchCallStack[pathSearchCallStackSize-2].nodeEdge)){
                        spqr_node head = findEdgeHeadNoCompression(dec, dfsData->nodeEdge);
                        spqr_node tail = findEdgeTailNoCompression(dec, dfsData->nodeEdge);
                        spqr_node other = head == dfsData->node ? tail : head;
                        assert(other != dfsData->node);
                        assert(!nodeVisited[other]);
                        if(other == target){
                            break;
                        }
                        //We go up a level: add new node to the call stack

                        pathSearchCallStack[pathSearchCallStackSize].node = other;
                        pathSearchCallStack[pathSearchCallStackSize].nodeEdge = getFirstNodeEdge(dec,other);
                        ++pathSearchCallStackSize;
                        continue;
                    }
                    do{
                        dfsData->nodeEdge = getNextNodeEdgeNoCompression(dec,dfsData->nodeEdge,dfsData->node);
                        if(dfsData->nodeEdge == getFirstNodeEdge(dec,dfsData->node)){
                            --pathSearchCallStackSize;
                            dfsData = &pathSearchCallStack[pathSearchCallStackSize-1];
                        }else{
                            break;
                        }
                    }while(pathSearchCallStackSize > 0);
                }
                for (int i = 0; i < pathSearchCallStackSize; ++i) {
                    if(edgeIsTree(dec,pathSearchCallStack[i].nodeEdge)){
                        process_edge(output,&num_rows,callStack,&callStackSize,pathSearchCallStack[i].nodeEdge,dec);
                    }
                }

                pathSearchCallStackSize = 0;
                break;
            }
            case SPQR_MEMBERTYPE_LOOP:
            case SPQR_MEMBERTYPE_PARALLEL:
            {
                spqr_edge first_edge = getFirstMemberEdge(dec, column_edge_member);
                spqr_edge iter_edge = first_edge;
                int tree_count = 0;
                do
                {
                    if(edgeIsTree(dec,iter_edge)){
                        process_edge(output,&num_rows,callStack,&callStackSize,iter_edge,dec);
                        tree_count++;
                    }
                    iter_edge = getNextMemberEdge(dec,iter_edge);
                }
                while(iter_edge != first_edge);
                if(tree_count > 1){
                    return -1;
                }
                break;
            }
            case SPQR_MEMBERTYPE_SERIES:
            {
                spqr_edge first_edge = getFirstMemberEdge(dec, column_edge_member);
                spqr_edge iter_edge = first_edge;
                int nontree_count = 0;
                do
                {
                    if(edgeIsTree(dec,iter_edge)){
                        process_edge(output,&num_rows,callStack,&callStackSize,iter_edge,dec);
                    }else{
                        nontree_count++;
                    }
                    iter_edge = getNextMemberEdge(dec,iter_edge);
                }
                while(iter_edge != first_edge);
                if(nontree_count != 1){
                    return -1;
                }
                break;
            }
            case SPQR_MEMBERTYPE_UNASSIGNED:
                assert(false);
        }
    }
    MATRECfreeBlockArray(dec->env,&pathSearchCallStack);
    MATRECfreeBlockArray(dec->env,&nodeVisited);
    MATRECfreeBlockArray(dec->env,&callStack);
    return num_rows;
}

static int qsort_integer_comparison (const void * a, const void * b)
{
    int *s1 = (int *)a;
    int *s2 = (int *)b;

    if(*s1 > *s2) {
        return 1;
    }
    else if(*s1 == *s2) {
        return 0;
    }
    else {
        return -1;
    }
}
bool MATRECGraphicDecompositionVerifyCycle(const MATRECGraphicDecomposition * dec, MATREC_col column, MATREC_row * column_rows,
                                           int num_rows, MATREC_row * computed_column_storage){
    int num_found_rows = decompositionGetFundamentalCycleRows(dec,column,computed_column_storage);

    if(num_found_rows != num_rows){
        return false;
    }
    if(num_rows == 0){
        return true;
    }
    qsort(computed_column_storage, (size_t) num_rows, sizeof(MATREC_row), qsort_integer_comparison);
    qsort(column_rows            , (size_t) num_rows, sizeof(MATREC_row), qsort_integer_comparison);

    for (int i = 0; i < num_rows; ++i) {
        if(column_rows[i] != computed_column_storage[i]){
            return false;
        }
    }
    return true;
}

static spqr_member largestMemberID(const MATRECGraphicDecomposition *dec){
    return dec->numMembers;
}
static spqr_edge largestEdgeID(const MATRECGraphicDecomposition *dec){
    return dec->numEdges;
}
static spqr_node largestNodeID(const MATRECGraphicDecomposition *dec){
    return dec->numNodes;
}
static int numConnectedComponents(const MATRECGraphicDecomposition *dec){
    return dec->numConnectedComponents;
}
static MATREC_ERROR createChildMarker(MATRECGraphicDecomposition *dec, spqr_member member, spqr_member child, bool isTree, spqr_edge * pEdge){
    MATREC_CALL(createEdge(dec,member,pEdge));
    dec->edges[*pEdge].element = isTree ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;
    dec->edges[*pEdge].childMember = child;

    addEdgeToMemberEdgeList(dec,*pEdge,member);
    return MATREC_OKAY;
}
static MATREC_ERROR createParentMarker(MATRECGraphicDecomposition *dec, spqr_member member, bool isTree, spqr_member parent, spqr_member parentMarker
        , spqr_edge * edge){

    MATREC_CALL(createEdge(dec,member,edge));
    dec->edges[*edge].element = isTree ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;

    addEdgeToMemberEdgeList(dec,*edge,member);

    dec->members[member].parentMember = parent;
    dec->members[member].markerOfParent = parentMarker;
    dec->members[member].markerToParent = *edge;
    return MATREC_OKAY;
}
static MATREC_ERROR createMarkerPair(MATRECGraphicDecomposition *dec, spqr_member parentMember, spqr_member childMember, bool parentIsTree){
    spqr_edge parentToChildMarker = SPQR_INVALID_EDGE;
    MATREC_CALL(createChildMarker(dec,parentMember,childMember,parentIsTree,&parentToChildMarker));

    spqr_edge childToParentMarker = SPQR_INVALID_EDGE;
    MATREC_CALL(createParentMarker(dec,childMember,!parentIsTree,parentMember,parentToChildMarker,&childToParentMarker));

    return MATREC_OKAY;
}
static MATREC_ERROR createMarkerPairWithReferences(MATRECGraphicDecomposition *dec, spqr_member parentMember, spqr_member childMember, bool parentIsTree,
                                          spqr_edge * parentToChild, spqr_edge *childToParent){
    MATREC_CALL(createChildMarker(dec,parentMember,childMember,parentIsTree,parentToChild));
    MATREC_CALL(createParentMarker(dec,childMember,!parentIsTree,parentMember,*parentToChild,childToParent));

    return MATREC_OKAY;
}

static void moveEdgeToNewMember(MATRECGraphicDecomposition *dec, spqr_edge edge, spqr_member oldMember, spqr_member newMember){
    assert(SPQRedgeIsValid(edge));
    assert(edge < dec->memEdges);
    assert(dec);

    assert(memberIsRepresentative(dec,oldMember));
    assert(memberIsRepresentative(dec,newMember));
    //Need to change the edge's member, remove it from the current member list and add it to the new member list
    assert(findEdgeMemberNoCompression(dec,edge) == oldMember);

    removeEdgeFromMemberEdgeList(dec,edge,oldMember);
    addEdgeToMemberEdgeList(dec,edge,newMember);

    dec->edges[edge].member = newMember;

    //If this edge has a childMember, update the information correctly!
    spqr_member childMember = dec->edges[edge].childMember;
    if(SPQRmemberIsValid(childMember)){
        spqr_member childRepresentative = findEdgeChildMember(dec, edge);
        dec->members[childRepresentative].parentMember = newMember;
    }
    //If this edge is a marker to the parent, update the child edge marker of the parent to reflect the move
    if(dec->members[oldMember].markerToParent == edge){
        dec->members[newMember].markerToParent = edge;
        dec->members[newMember].parentMember = dec->members[oldMember].parentMember;
        dec->members[newMember].markerOfParent = dec->members[oldMember].markerOfParent;

        assert(findEdgeChildMemberNoCompression(dec,dec->members[oldMember].markerOfParent) == oldMember);
        dec->edges[dec->members[oldMember].markerOfParent].childMember = newMember;
    }
}
static void mergeMemberEdgeList(MATRECGraphicDecomposition *dec, spqr_member toMergeInto, spqr_member toRemove){
    spqr_edge firstIntoEdge = getFirstMemberEdge(dec, toMergeInto);
    spqr_edge firstFromEdge = getFirstMemberEdge(dec, toRemove);
    assert(SPQRedgeIsValid(firstIntoEdge));
    assert(SPQRedgeIsValid(firstFromEdge));

    spqr_edge lastIntoEdge = getPreviousMemberEdge(dec, firstIntoEdge);
    spqr_edge lastFromEdge = getPreviousMemberEdge(dec, firstFromEdge);

    //Relink linked lists to merge them effectively
    dec->edges[firstIntoEdge].edgeListNode.previous = lastFromEdge;
    dec->edges[lastIntoEdge].edgeListNode.next = firstFromEdge;
    dec->edges[firstFromEdge].edgeListNode.previous = lastIntoEdge;
    dec->edges[lastFromEdge].edgeListNode.next = firstIntoEdge;

    //Clean up old
    dec->members[toMergeInto].num_edges += dec->members[toRemove].num_edges;
    dec->members[toRemove].num_edges = 0;
    dec->members[toRemove].firstEdge = SPQR_INVALID_EDGE;

}

static void changeLoopToSeries(MATRECGraphicDecomposition * dec, spqr_member member){
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(dec);
    assert((getMemberType(dec,member) == SPQR_MEMBERTYPE_PARALLEL || getMemberType(dec, member) == SPQR_MEMBERTYPE_SERIES || getMemberType(dec,member) == SPQR_MEMBERTYPE_LOOP) && getNumMemberEdges(dec, member) == 2);
    assert(memberIsRepresentative(dec,member));
    dec->members[member].type = SPQR_MEMBERTYPE_SERIES;
}
static void changeLoopToParallel(MATRECGraphicDecomposition * dec, spqr_member member){
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(dec);
    assert((getMemberType(dec,member) == SPQR_MEMBERTYPE_PARALLEL
    || getMemberType(dec, member) == SPQR_MEMBERTYPE_SERIES
    || getMemberType(dec,member) == SPQR_MEMBERTYPE_LOOP) && getNumMemberEdges(dec, member) == 2);
    assert(memberIsRepresentative(dec,member));
    dec->members[member].type = SPQR_MEMBERTYPE_PARALLEL;
}
bool MATRECGraphicDecompositionIsMinimal(const MATRECGraphicDecomposition * dec){
    //Relies on parents/children etc. being set correctly in the tree
    bool isMinimal = true;
    for (spqr_member member = 0; member < dec->numMembers; ++member) {
        if (!memberIsRepresentative(dec, member) || getMemberType(dec,member) == SPQR_MEMBERTYPE_UNASSIGNED ){
            continue;
        } //TODO: fix loop making INVALID members here... this is not a pretty way to do this
        spqr_member memberParent = findMemberParentNoCompression(dec, member);
        if(SPQRmemberIsValid(memberParent)){
            SPQRMemberType memberType = getMemberType(dec, member);
            SPQRMemberType parentType = getMemberType(dec, memberParent);
            if(memberType == parentType && memberType != SPQR_MEMBERTYPE_RIGID){
                isMinimal = false;
                break;
            }
        }

    }
    return isMinimal;
}

static void mergeLoop(MATRECGraphicDecomposition * dec, spqr_member loopMember){
    assert(getMemberType(dec,loopMember) == SPQR_MEMBERTYPE_SERIES ||
    getMemberType(dec, loopMember) == SPQR_MEMBERTYPE_PARALLEL);
    assert(getNumMemberEdges(dec,loopMember) == 2);
    assert(memberIsRepresentative(dec,loopMember));

#ifndef NDEBUG
    {
        int num_markers = 0;
        spqr_edge first_edge = getFirstMemberEdge(dec, loopMember);
        spqr_edge edge = first_edge;
        do {
            if(edgeIsMarker(dec,edge) || (markerToParent(dec,loopMember) == edge)){
                num_markers++;
            }
            edge = getNextMemberEdge(dec, edge);
        } while (edge != first_edge);
        assert(num_markers == 2);
    };
#endif

    spqr_member seriesMember = SPQR_INVALID_MEMBER;
    spqr_edge loopSeriesEdge = SPQR_INVALID_EDGE;
    bool seriesIsParent = false;
    spqr_member otherMember = SPQR_INVALID_MEMBER;
    spqr_edge loopOtherEdge = SPQR_INVALID_EDGE;
    bool otherIsParent = false;
    {
        spqr_edge first_edge = getFirstMemberEdge(dec, loopMember);
        spqr_edge edge = first_edge;
        do {
            if (edgeIsMarker(dec, edge)) {
                spqr_member child = findEdgeChildMember(dec, edge);
                if(getMemberType(dec,child) == SPQR_MEMBERTYPE_SERIES){
                    seriesMember = child;
                    loopSeriesEdge = edge;
                    seriesIsParent = false;
                }else{
                    otherMember = child;
                    loopOtherEdge = edge;
                    otherIsParent = false;
                }
            }else if (markerToParent(dec,loopMember) == edge){
                spqr_member parent = findMemberParent(dec,loopMember);
                if(getMemberType(dec,parent) == SPQR_MEMBERTYPE_SERIES){
                    seriesMember = parent;
                    loopSeriesEdge = edge;
                    seriesIsParent = true;
                }else{
                    otherMember = parent;
                    loopOtherEdge = edge;
                    otherIsParent = true;
                }
            }
            edge = getNextMemberEdge(dec, edge);
        } while (edge != first_edge);
    };
    assert(SPQRmemberIsValid(seriesMember) && SPQRmemberIsValid(otherMember));
    assert(getMemberType(dec,seriesMember) == SPQR_MEMBERTYPE_SERIES &&
            getMemberType(dec,otherMember) != SPQR_MEMBERTYPE_SERIES);
    assert(!(seriesIsParent && otherIsParent));

    if(seriesIsParent){
        //other member must be a child
        spqr_member seriesChildEdge = markerOfParent(dec,loopMember);
        dec->members[otherMember].markerOfParent = seriesChildEdge;
        dec->members[otherMember].parentMember = seriesMember;
        dec->edges[seriesChildEdge].childMember = otherMember;

        removeEdgeFromMemberEdgeList(dec,loopSeriesEdge,loopMember);
        removeEdgeFromMemberEdgeList(dec,loopOtherEdge,loopMember);
        dec->members[loopMember].type = SPQR_MEMBERTYPE_UNASSIGNED;
    }else if(otherIsParent){
        //series member is a child
        spqr_member otherChildEdge = markerOfParent(dec,loopMember);
        dec->members[seriesMember].markerOfParent = otherChildEdge;
        dec->members[seriesMember].parentMember = otherMember;
        dec->edges[otherChildEdge].childMember = seriesMember;

        removeEdgeFromMemberEdgeList(dec,loopSeriesEdge,loopMember);
        removeEdgeFromMemberEdgeList(dec,loopOtherEdge,loopMember);
        dec->members[loopMember].type = SPQR_MEMBERTYPE_UNASSIGNED;
    }else{
        //The loop member is the root; we make the new series member the root
        spqr_edge seriesArcToLoop = markerToParent(dec,seriesMember);

        dec->members[seriesMember].markerOfParent = SPQR_INVALID_EDGE;
        dec->members[seriesMember].markerToParent = SPQR_INVALID_EDGE;
        dec->members[seriesMember].parentMember = SPQR_INVALID_MEMBER;
        dec->edges[seriesArcToLoop].childMember = otherMember;

        dec->members[otherMember].parentMember = seriesMember;
        dec->members[otherMember].markerOfParent = seriesArcToLoop;

        removeEdgeFromMemberEdgeList(dec,loopSeriesEdge,loopMember);
        removeEdgeFromMemberEdgeList(dec,loopOtherEdge,loopMember);
        dec->members[loopMember].type = SPQR_MEMBERTYPE_UNASSIGNED;
    }


}

static void decreaseNumConnectedComponents(MATRECGraphicDecomposition *dec, int by){
    dec->numConnectedComponents-= by;
    assert(dec->numConnectedComponents >= 1);
}

static void reorderComponent(MATRECGraphicDecomposition *dec, spqr_member newRoot){
    assert(dec);
    assert(memberIsRepresentative(dec,newRoot));
    //If the newRoot has no parent, it is already the root, so then there's no need to reorder.
    if(SPQRmemberIsValid(dec->members[newRoot].parentMember)){
        spqr_member member = findMemberParent(dec, newRoot);
        spqr_member newParent = newRoot;
        spqr_edge newMarkerToParent = dec->members[newRoot].markerOfParent;
        spqr_edge markerOfNewParent = dec->members[newRoot].markerToParent;

        //Recursively update the parent
        do{
            assert(SPQRmemberIsValid(member));
            assert(SPQRmemberIsValid(newParent));
            spqr_member oldParent = findMemberParent(dec, member);
            spqr_edge oldMarkerToParent = dec->members[member].markerToParent;
            spqr_edge oldMarkerOfParent = dec->members[member].markerOfParent;

            dec->members[member].markerToParent = newMarkerToParent;
            dec->members[member].markerOfParent = markerOfNewParent;
            dec->members[member].parentMember = newParent;
            dec->edges[markerOfNewParent].childMember = member;
            dec->edges[newMarkerToParent].childMember = -1;

            if (SPQRmemberIsValid(oldParent)){
                newParent = member;
                member = oldParent;
                newMarkerToParent = oldMarkerOfParent;
                markerOfNewParent = oldMarkerToParent;
            }else{
                break;
            }
        }while(true);
        dec->members[newRoot].parentMember = SPQR_INVALID_MEMBER;
        dec->members[newRoot].markerToParent = SPQR_INVALID_EDGE;
        dec->members[newRoot].markerOfParent = SPQR_INVALID_EDGE;
    }
}

static void edgeToDot(FILE * stream, const MATRECGraphicDecomposition * dec,
                      spqr_edge edge, unsigned long dot_head, unsigned long dot_tail, bool useElementNames){
    assert(SPQRedgeIsValid(edge));
    spqr_member member = findEdgeMemberNoCompression(dec, edge);
    SPQRMemberType member_type = getMemberType(dec, member);
    char type = typeToChar(member_type);
    const char* color = edgeIsTree(dec,edge) ? ",color=red" :",color=blue";

    int edge_name = edge;

    if(markerToParent(dec,member) == edge){
        if(useElementNames){
            edge_name = -1;
        }
        fprintf(stream, "    %c_%d_%lu -> %c_p_%d [label=\"%d\",style=dashed%s];\n", type, member, dot_head, type, member, edge_name, color);
        fprintf(stream, "    %c_p_%d -> %c_%d_%lu [label=\"%d\",style=dashed%s];\n", type, member, type, member, dot_tail, edge_name, color);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
        fprintf(stream, "    %c_p_%d [style=dashed];\n", type, member);
    }else if(edgeIsMarker(dec,edge)){
        spqr_member child = findEdgeChildMemberNoCompression(dec, edge);
        char childType = typeToChar(getMemberType(dec,child));
        if(useElementNames){
            edge_name = -1;
        }
        fprintf(stream, "    %c_%d_%lu -> %c_c_%d [label=\"%d\",style=dotted%s];\n", type, member, dot_head, type, child, edge_name, color);
        fprintf(stream, "    %c_c_%d -> %c_%d_%lu [label=\"%d\",style=dotted%s];\n", type, child, type, member, dot_tail, edge_name, color);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
        fprintf(stream, "    %c_c_%d [style=dotted];\n", type, child);
        fprintf(stream, "    %c_p_%d -> %c_c_%d [style=dashed,dir=forward];\n", childType, child, type, child);
    }else{
        if(useElementNames){
            spqr_element element = dec->edges[edge].element;
            if(SPQRelementIsRow(element)){
                edge_name = (int) SPQRelementToRow(element);
            }else{
                edge_name = (int) SPQRelementToColumn(element);
            }
        }

        fprintf(stream, "    %c_%d_%lu -> %c_%d_%lu [label=\"%d \",style=bold%s];\n", type, member, dot_head, type, member, dot_tail,
                edge_name, color);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
    }
}

static void decompositionToDot(FILE * stream, const MATRECGraphicDecomposition *dec, bool useElementNames ){
    fprintf(stream, "//decomposition\ndigraph decomposition{\n   compound = true;\n");
    for (spqr_member member = 0; member < dec->numMembers; ++member){
        if(!memberIsRepresentative(dec,member)) continue;
        fprintf(stream,"   subgraph member_%d{\n",member);
        switch(getMemberType(dec,member)){
            case SPQR_MEMBERTYPE_RIGID:
            {
                spqr_edge first_edge = getFirstMemberEdge(dec, member);
                spqr_edge edge = first_edge;
                do
                {
                    unsigned long edgeHead = (unsigned long) findEdgeHeadNoCompression(dec,edge);
                    unsigned long edgeTail = (unsigned long) findEdgeTailNoCompression(dec,edge);
                    edgeToDot(stream,dec,edge,edgeHead,edgeTail,useElementNames);
                    edge = getNextMemberEdge(dec,edge);
                }
                while(edge != first_edge);
                break;
            }
            case SPQR_MEMBERTYPE_LOOP:
            case SPQR_MEMBERTYPE_PARALLEL:
            {
                spqr_edge first_edge = getFirstMemberEdge(dec, member);
                spqr_edge edge = first_edge;
                do
                {
                    edgeToDot(stream,dec,edge,0,1,useElementNames);
                    edge = getNextMemberEdge(dec,edge);
                }
                while(edge != first_edge);
                break;
            }
            case SPQR_MEMBERTYPE_SERIES:
            {
                unsigned long i = 0;
                unsigned long num_member_edges = (unsigned long) getNumMemberEdges(dec, member);
                spqr_edge first_edge = getFirstMemberEdge(dec, member);
                spqr_edge edge = first_edge;
                do {
                    edgeToDot(stream, dec, edge, i, (i + 1) % num_member_edges,useElementNames);
                    edge = getNextMemberEdge(dec, edge);
                    i++;
                } while (edge != first_edge);
                break;
            }
            case SPQR_MEMBERTYPE_UNASSIGNED:
                break;
        }
        fprintf(stream,"   }\n");
    }
    fprintf(stream,"}\n");
}

static int max(int a, int b){
    return (a > b) ? a : b;
}

typedef int path_edge_id;
#define INVALID_PATH_EDGE (-1)

static bool pathEdgeIsInvalid(const path_edge_id edge) {
    return edge < 0;
}

static bool pathEdgeIsValid(const path_edge_id edge) {
    return !pathEdgeIsInvalid(edge);
}

typedef struct {
    spqr_edge edge;
    spqr_node edgeHead; //These can be used in various places to prevent additional find()'s
    spqr_node edgeTail;
    path_edge_id nextMember;
    path_edge_id nextOverall;
} PathEdgeListNode;

typedef int reduced_member_id;
#define INVALID_REDUCED_MEMBER (-1)

static bool reducedMemberIsInvalid(const reduced_member_id id) {
    return id < 0;
}
static bool reducedMemberIsValid(const reduced_member_id id){
    return !reducedMemberIsInvalid(id);
}

typedef int children_idx;

typedef enum {
    TYPE_INVALID = 1,
    TYPE_SINGLE_CHILD = 2,
    TYPE_DOUBLE_CHILD = 3,
    TYPE_CYCLE_CHILD = 4,
    TYPE_ROOT = 5
    //TODO fix
} ReducedMemberType;

typedef struct {
    spqr_member member;
    spqr_member rootMember;
    int depth;
    ReducedMemberType type;
    reduced_member_id parent;

    children_idx firstChild;
    children_idx numChildren;

    path_edge_id firstPathEdge;
    int numPathEdges;

    int numOneEnd;
    int numTwoEnds;
    spqr_edge childMarkerEdges[2];
    spqr_node rigidEndNodes[4];
} MATRECColReducedMember;

typedef struct {
    int rootDepth;
    reduced_member_id root;
} MATRECColReducedComponent;

typedef struct {
    reduced_member_id reducedMember;
    reduced_member_id rootDepthMinimizer;
} MemberInfo;

typedef struct {
    spqr_member member;
} CreateReducedMembersCallstack;

struct MATRECGraphicColumnAdditionImpl {
    bool remainsGraphic;

    MATRECColReducedMember *reducedMembers;
    int memReducedMembers;
    int numReducedMembers;

    MATRECColReducedComponent *reducedComponents;
    int memReducedComponents;
    int numReducedComponents;

    MemberInfo *memberInformation;
    int memMemberInformation;
    int numMemberInformation;

    reduced_member_id *childrenStorage;
    int memChildrenStorage;
    int numChildrenStorage;

    PathEdgeListNode *pathEdges;
    int memPathEdges;
    int numPathEdges;
    path_edge_id firstOverallPathEdge;

    int *nodePathDegree;
    int memNodePathDegree;

    bool *edgeInPath;
    int memEdgesInPath;

    CreateReducedMembersCallstack * createReducedMembersCallStack;
    int memCreateReducedMembersCallStack;

    MATREC_col newColIndex;

    MATREC_row *newRowEdges;
    int memNewRowEdges;
    int numNewRowEdges;

    spqr_edge *decompositionRowEdges;
    int memDecompositionRowEdges;
    int numDecompositionRowEdges;
};

static void cleanupPreviousIteration(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol) {
    assert(dec);
    assert(newCol);

    path_edge_id pathEdge = newCol->firstOverallPathEdge;
    while (pathEdgeIsValid(pathEdge)) {
        spqr_node head = newCol->pathEdges[pathEdge].edgeHead;
        spqr_node tail = newCol->pathEdges[pathEdge].edgeTail;
        if(SPQRnodeIsValid(head)){
            newCol->nodePathDegree[head] = 0;
        }
        if(SPQRnodeIsValid(tail)){
            newCol->nodePathDegree[tail] = 0;
        }

        spqr_edge edge = newCol->pathEdges[pathEdge].edge;
        if(edge < newCol->memEdgesInPath){
            newCol->edgeInPath[edge] = false;
        }
        pathEdge = newCol->pathEdges[pathEdge].nextOverall;
    }
#ifndef NDEBUG
    for (int i = 0; i < newCol->memEdgesInPath; ++i) {
        assert(newCol->edgeInPath[i] == false);
    }

    for (int i = 0; i < newCol->memNodePathDegree; ++i) {
        assert(newCol->nodePathDegree[i] == 0);
    }
#endif

    newCol->firstOverallPathEdge = INVALID_PATH_EDGE;
    newCol->numPathEdges = 0;
}

MATREC_ERROR MATRECcreateGraphicColumnAddition(MATREC *env, MATRECGraphicColumnAddition **pNewCol) {
    assert(env);

    MATREC_CALL(MATRECallocBlock(env, pNewCol));
    MATRECGraphicColumnAddition *newCol = *pNewCol;

    newCol->remainsGraphic = false;
    newCol->reducedMembers = NULL;
    newCol->memReducedMembers = 0;
    newCol->numReducedMembers = 0;

    newCol->reducedComponents = NULL;
    newCol->memReducedComponents = 0;
    newCol->numReducedComponents = 0;

    newCol->memberInformation = NULL;
    newCol->memMemberInformation = 0;
    newCol->numMemberInformation = 0;

    newCol->childrenStorage = NULL;
    newCol->memChildrenStorage = 0;
    newCol->numChildrenStorage = 0;

    newCol->pathEdges = NULL;
    newCol->memPathEdges = 0;
    newCol->numPathEdges = 0;
    newCol->firstOverallPathEdge = INVALID_PATH_EDGE;

    newCol->nodePathDegree = NULL;
    newCol->memNodePathDegree = 0;

    newCol->edgeInPath = NULL;
    newCol->memEdgesInPath = 0;

    newCol->createReducedMembersCallStack = NULL;
    newCol->memCreateReducedMembersCallStack = 0;

    newCol->newColIndex = MATREC_INVALID_COL;

    newCol->newRowEdges = NULL;
    newCol->memNewRowEdges = 0;
    newCol->numNewRowEdges = 0;

    newCol->decompositionRowEdges = NULL;
    newCol->memDecompositionRowEdges = 0;
    newCol->numDecompositionRowEdges = 0;

    return MATREC_OKAY;
}

void MATRECfreeGraphicColumnAddition(MATREC *env, MATRECGraphicColumnAddition **pNewCol) {
    assert(env);
    MATRECGraphicColumnAddition *newCol = *pNewCol;
    MATRECfreeBlockArray(env, &newCol->decompositionRowEdges);
    MATRECfreeBlockArray(env, &newCol->newRowEdges);
    MATRECfreeBlockArray(env, &newCol->createReducedMembersCallStack);
    MATRECfreeBlockArray(env, &newCol->edgeInPath);
    MATRECfreeBlockArray(env, &newCol->nodePathDegree);
    MATRECfreeBlockArray(env, &newCol->pathEdges);
    MATRECfreeBlockArray(env, &newCol->childrenStorage);
    MATRECfreeBlockArray(env, &newCol->memberInformation);
    MATRECfreeBlockArray(env, &newCol->reducedComponents);
    MATRECfreeBlockArray(env, &newCol->reducedMembers);

    MATRECfreeBlock(env, pNewCol);
}


static reduced_member_id createReducedMembersToRoot(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition * newCol, const spqr_member firstMember ){
    assert(SPQRmemberIsValid(firstMember));

    CreateReducedMembersCallstack * callstack = newCol->createReducedMembersCallStack;
    callstack[0].member = firstMember;
    int callDepth = 0;

    while(callDepth >= 0){
        spqr_member member = callstack[callDepth].member;
        reduced_member_id reducedMember = newCol->memberInformation[member].reducedMember;

        bool reducedValid = reducedMemberIsValid(reducedMember);
        if(!reducedValid) {
            //reduced member was not yet created; we create it
            reducedMember = newCol->numReducedMembers;

            MATRECColReducedMember *reducedMemberData = &newCol->reducedMembers[reducedMember];
            ++newCol->numReducedMembers;

            reducedMemberData->member = member;
            reducedMemberData->numChildren = 0;

            reducedMemberData->type = TYPE_INVALID;
            reducedMemberData->firstPathEdge = INVALID_PATH_EDGE;
            reducedMemberData->numPathEdges = 0;
            for (int i = 0; i < 4; ++i) {
                reducedMemberData->rigidEndNodes[i] = SPQR_INVALID_NODE;
            }
            //The children are set later

            newCol->memberInformation[member].reducedMember = reducedMember;
            assert(memberIsRepresentative(dec, member));
            spqr_member parentMember = findMemberParent(dec, member);

            if (SPQRmemberIsValid(parentMember)) {
                //recursive call to parent member
                ++callDepth;
                assert(callDepth < newCol->memCreateReducedMembersCallStack);
                callstack[callDepth].member = parentMember;
                continue;

            } else {
                //we found a new reduced decomposition component

                reducedMemberData->parent = INVALID_REDUCED_MEMBER;
                reducedMemberData->depth = 0;
                reducedMemberData->rootMember = member;

                assert(newCol->numReducedComponents < newCol->memReducedComponents);
                newCol->reducedComponents[newCol->numReducedComponents].root = reducedMember;
                ++newCol->numReducedComponents;
            }
        }
        if(reducedValid){
            assert(reducedMember < newCol->numReducedMembers);
            //Reduced member was already created in earlier call
            //update the depth of the root if appropriate
            reduced_member_id * depthMinimizer = &newCol->memberInformation[newCol->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
            if(reducedMemberIsInvalid(*depthMinimizer) ||
               newCol->reducedMembers[reducedMember].depth < newCol->reducedMembers[*depthMinimizer].depth){
                *depthMinimizer = reducedMember;
            }
        }
        while(true){
            --callDepth;
            if(callDepth < 0 ) break;
            spqr_member parentMember = callstack[callDepth + 1].member;
            reduced_member_id parentReducedMember = newCol->memberInformation[parentMember].reducedMember;
            spqr_member currentMember = callstack[callDepth].member;
            reduced_member_id currentReducedMember = newCol->memberInformation[currentMember].reducedMember;

            MATRECColReducedMember *parentReducedMemberData = &newCol->reducedMembers[parentReducedMember];
            MATRECColReducedMember *reducedMemberData = &newCol->reducedMembers[currentReducedMember];

            reducedMemberData->parent = parentReducedMember;
            reducedMemberData->depth = parentReducedMemberData->depth + 1;
            reducedMemberData->rootMember = parentReducedMemberData->rootMember;

            newCol->reducedMembers[parentReducedMember].numChildren++;
        }

    }

    reduced_member_id returnedMember = newCol->memberInformation[callstack[0].member].reducedMember;
    return returnedMember;
}

static MATREC_ERROR constructReducedDecomposition(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol) {
    assert(dec);
    assert(newCol);
#ifndef NDEBUG
    for (int i = 0; i < newCol->memMemberInformation; ++i) {
        assert(reducedMemberIsInvalid(newCol->memberInformation[i].reducedMember));
    }
#endif
    newCol->numReducedComponents = 0;
    newCol->numReducedMembers = 0;
    if (newCol->numDecompositionRowEdges == 0) { //Early return in case the reduced decomposition will be empty
        return MATREC_OKAY;
    }
    assert(newCol->numReducedMembers == 0);
    assert(newCol->numReducedComponents == 0);

    int newSize = largestMemberID(dec); //Is this sufficient?
    if (newSize > newCol->memReducedMembers) {
        newCol->memReducedMembers = max(2 * newCol->memReducedMembers, newSize);
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->reducedMembers, (size_t) newCol->memReducedMembers));
    }
    if (newSize > newCol->memMemberInformation) {
        int updatedSize = max(2 * newCol->memMemberInformation, newSize);
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->memberInformation, (size_t) updatedSize));
        for (int i = newCol->memMemberInformation; i < updatedSize; ++i) {
            newCol->memberInformation[i].reducedMember = INVALID_REDUCED_MEMBER;
            newCol->memberInformation[i].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
        }
        newCol->memMemberInformation = updatedSize;

    }

    int numComponents = numConnectedComponents(dec);
    if (numComponents > newCol->memReducedComponents) {
        newCol->memReducedComponents = max(2 * newCol->memReducedComponents, numComponents);
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->reducedComponents, (size_t) newCol->memReducedComponents));
    }

    int numMembers = getNumMembers(dec);
    if (newCol->memCreateReducedMembersCallStack < numMembers) {
        newCol->memCreateReducedMembersCallStack = max(2 * newCol->memCreateReducedMembersCallStack, numMembers);
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->createReducedMembersCallStack,
                                        (size_t) newCol->memCreateReducedMembersCallStack));
    }

    //Create the reduced members (recursively)
    for (int i = 0; i < newCol->numDecompositionRowEdges; ++i) {
        assert(i < newCol->memDecompositionRowEdges);
        spqr_edge edge = newCol->decompositionRowEdges[i];
        spqr_member edgeMember = findEdgeMember(dec, edge);
        reduced_member_id reducedMember = createReducedMembersToRoot(dec, newCol, edgeMember);
        reduced_member_id *depthMinimizer = &newCol->memberInformation[newCol->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
        if (reducedMemberIsInvalid(*depthMinimizer)) {
            *depthMinimizer = reducedMember;
        }
    }

    //Set the reduced roots according to the root depth minimizers
    for (int i = 0; i < newCol->numReducedComponents; ++i) {
        MATRECColReducedComponent *component = &newCol->reducedComponents[i];
        spqr_member rootMember = newCol->reducedMembers[component->root].member;
        reduced_member_id reducedMinimizer = newCol->memberInformation[rootMember].rootDepthMinimizer;
        component->rootDepth = newCol->reducedMembers[reducedMinimizer].depth;
        component->root = reducedMinimizer;

        //This simplifies code further down which does not need to be component-aware; just pretend that the reduced member is the new root.
        newCol->reducedMembers[component->root].parent = INVALID_REDUCED_MEMBER;
        assert(memberIsRepresentative(dec, rootMember));
    }

    //update the children array
    int numTotalChildren = 0;
    for (int i = 0; i < newCol->numReducedMembers; ++i) {
        MATRECColReducedMember *reducedMember = &newCol->reducedMembers[i];
        reduced_member_id minimizer = newCol->memberInformation[reducedMember->rootMember].rootDepthMinimizer;
        if (reducedMember->depth >= newCol->reducedMembers[minimizer].depth) {
            reducedMember->firstChild = numTotalChildren;
            numTotalChildren += reducedMember->numChildren;
            reducedMember->numChildren = 0;
        }
    }

    if (newCol->memChildrenStorage < numTotalChildren) {
        int newMemSize = max(newCol->memChildrenStorage * 2, numTotalChildren);
        newCol->memChildrenStorage = newMemSize;
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->childrenStorage, (size_t) newCol->memChildrenStorage));
    }
    newCol->numChildrenStorage = numTotalChildren;

    //Fill up the children array`
    for (reduced_member_id reducedMember = 0; reducedMember < newCol->numReducedMembers; ++reducedMember) {
        MATRECColReducedMember *reducedMemberData = &newCol->reducedMembers[reducedMember];
        if (reducedMemberData->depth <=
            newCol->reducedMembers[newCol->memberInformation[reducedMemberData->rootMember].rootDepthMinimizer].depth) {
            continue;
        }
        spqr_member parentMember = findMemberParent(dec, reducedMemberData->member);
        reduced_member_id parentReducedMember = SPQRmemberIsValid(parentMember)
                                                ? newCol->memberInformation[parentMember].reducedMember
                                                : INVALID_REDUCED_MEMBER;
        if (reducedMemberIsValid(parentReducedMember)) {
            //TODO: probably one of these two checks/branches is unnecessary, as there is a single failure case? (Not sure)
            MATRECColReducedMember *parentReducedMemberData = &newCol->reducedMembers[parentReducedMember];
            newCol->childrenStorage[parentReducedMemberData->firstChild +
                                    parentReducedMemberData->numChildren] = reducedMember;
            ++parentReducedMemberData->numChildren;
        }
    }

    //Clean up the root depth minimizers.
    for (int i = 0; i < newCol->numReducedMembers; ++i) {
        MATRECColReducedMember *reducedMember = &newCol->reducedMembers[i];
        assert(reducedMember);
        spqr_member rootMember = reducedMember->rootMember;
        assert(rootMember >= 0);
        assert(rootMember < dec->memMembers);
        newCol->memberInformation[rootMember].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
    }

    return MATREC_OKAY;
}

static void cleanUpMemberInformation(MATRECGraphicColumnAddition * newCol){
    //This loop is at the end as memberInformation is also used to assign the cut edges during propagation
    //Clean up the memberInformation array
    for (int i = 0; i < newCol->numReducedMembers; ++i) {
        newCol->memberInformation[newCol->reducedMembers[i].member].reducedMember = INVALID_REDUCED_MEMBER;
    }
#ifndef NDEBUG
    for (int i = 0; i < newCol->memMemberInformation; ++i) {
        assert(reducedMemberIsInvalid(newCol->memberInformation[i].reducedMember));
    }
#endif
}

static void createPathEdge(
        MATRECGraphicDecomposition * dec, MATRECGraphicColumnAddition * newCol,
        const spqr_edge edge, const reduced_member_id reducedMember){
    assert(dec);
    assert(newCol);

    path_edge_id path_edge = newCol->numPathEdges;
    PathEdgeListNode * listNode = &newCol->pathEdges[path_edge];
    listNode->edge = edge;

    listNode->nextMember = newCol->reducedMembers[reducedMember].firstPathEdge;
    newCol->reducedMembers[reducedMember].firstPathEdge = path_edge;
    newCol->reducedMembers[reducedMember].numPathEdges += 1;

    listNode->nextOverall = newCol->firstOverallPathEdge;
    newCol->firstOverallPathEdge = path_edge;

    ++newCol->numPathEdges;
    assert(newCol->numPathEdges <= newCol->memPathEdges);

    assert(edge < newCol->memEdgesInPath);
    newCol->edgeInPath[edge] = true;
    if(getMemberType(dec,newCol->reducedMembers[reducedMember].member) == SPQR_MEMBERTYPE_RIGID){

        listNode->edgeHead = findEdgeHead(dec,edge);
        listNode->edgeTail = findEdgeTail(dec,edge);
        assert(SPQRnodeIsValid(listNode->edgeHead) && SPQRnodeIsValid(listNode->edgeTail));
        assert(listNode->edgeHead < newCol->memNodePathDegree && listNode->edgeTail < newCol->memNodePathDegree);
        ++newCol->nodePathDegree[listNode->edgeHead];
        ++newCol->nodePathDegree[listNode->edgeTail];
    }else{
        listNode->edgeHead = SPQR_INVALID_NODE;
        listNode->edgeTail = SPQR_INVALID_NODE;
    }

}

static MATREC_ERROR createPathEdges(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol){
    int maxNumPathEdges = newCol->numDecompositionRowEdges + getNumMembers(dec);
    if(newCol->memPathEdges < maxNumPathEdges){
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->pathEdges,(size_t) maxNumPathEdges)); //TODO: fix reallocation strategy
        newCol->memPathEdges = maxNumPathEdges;
    }
    int maxPathEdgeIndex = largestEdgeID(dec);
    if(newCol->memEdgesInPath < maxPathEdgeIndex){
        int newSize = maxPathEdgeIndex;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->edgeInPath,(size_t) newSize));//TODO: fix reallocation strategy
        for (int i = newCol->memEdgesInPath; i < newSize; ++i) {
            newCol->edgeInPath[i] = false;
        }
        newCol->memEdgesInPath = newSize;
    }
    int maxNumNodes = largestNodeID(dec);
    if(newCol->memNodePathDegree < maxNumNodes){
        int newSize = maxNumNodes;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->nodePathDegree,(size_t) newSize));
        for (int i = newCol->memNodePathDegree; i < newSize; ++i) {
            newCol->nodePathDegree[i] = 0;
        }
        newCol->memNodePathDegree = newSize;
    }
    for (int i = 0; i < newCol->numDecompositionRowEdges; ++i) {
        spqr_edge edge = newCol->decompositionRowEdges[i];
        spqr_member member = findEdgeMember(dec, edge);
        reduced_member_id reducedMember = newCol->memberInformation[member].reducedMember;
        createPathEdge(dec,newCol,edge,reducedMember);
    }

    return MATREC_OKAY;
}


/**
 * Saves the information of the current row and partitions it based on whether or not the given columns are
 * already part of the decomposition.
 */
static MATREC_ERROR
newColUpdateColInformation(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol, MATREC_col column, const MATREC_row *rows,
                           size_t numRows) {
    newCol->newColIndex = column;

    newCol->numDecompositionRowEdges = 0;
    newCol->numNewRowEdges = 0;

    for (size_t i = 0; i < numRows; ++i) {
        spqr_edge rowEdge = getDecompositionRowEdge(dec, rows[i]);
        if (SPQRedgeIsValid(rowEdge)) { //If the edge is the current decomposition: save it in the array
            if (newCol->numDecompositionRowEdges == newCol->memDecompositionRowEdges) {
                int newNumEdges = newCol->memDecompositionRowEdges == 0 ? 8 : 2 *
                                                                              newCol->memDecompositionRowEdges; //TODO: make reallocation numbers more consistent with rest?
                newCol->memDecompositionRowEdges = newNumEdges;
                MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->decompositionRowEdges,
                                                (size_t) newCol->memDecompositionRowEdges));
            }
            newCol->decompositionRowEdges[newCol->numDecompositionRowEdges] = rowEdge;
            ++newCol->numDecompositionRowEdges;
        } else {
            //Not in the decomposition: add it to the set of edges which are newly added with this row.
            if (newCol->numNewRowEdges == newCol->memNewRowEdges) {
                int newNumEdges = newCol->memNewRowEdges == 0 ? 8 : 2 *
                                                                    newCol->memNewRowEdges; //TODO: make reallocation numbers more consistent with rest?
                newCol->memNewRowEdges = newNumEdges;
                MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->newRowEdges,
                                                (size_t) newCol->memNewRowEdges));
            }
            newCol->newRowEdges[newCol->numNewRowEdges] = rows[i];
            newCol->numNewRowEdges++;
        }
    }

    return MATREC_OKAY;
}

static void countChildrenTypes(MATRECGraphicDecomposition* dec, MATRECGraphicColumnAddition * newCol, reduced_member_id reducedMember){
    newCol->reducedMembers[reducedMember].numOneEnd = 0;
    newCol->reducedMembers[reducedMember].numTwoEnds = 0;
    newCol->reducedMembers[reducedMember].childMarkerEdges[0] = SPQR_INVALID_EDGE;
    newCol->reducedMembers[reducedMember].childMarkerEdges[1] = SPQR_INVALID_EDGE;

    int nextChildMarker = 0;
    for (children_idx idx = newCol->reducedMembers[reducedMember].firstChild;
         idx < newCol->reducedMembers[reducedMember].firstChild
               + newCol->reducedMembers[reducedMember].numChildren;
         ++idx) {
        reduced_member_id reducedChild = newCol->childrenStorage[idx];
        assert(reducedMemberIsValid(reducedChild));
        if(newCol->reducedMembers[reducedChild].type == TYPE_SINGLE_CHILD){
            if(nextChildMarker < 2){
                newCol->reducedMembers[reducedMember].childMarkerEdges[nextChildMarker] = markerOfParent(dec, findMember(dec,newCol->reducedMembers[reducedChild].member)); //TODO: check if find is necessary
                ++nextChildMarker;
            }
            newCol->reducedMembers[reducedMember].numOneEnd++;
        }else if(newCol->reducedMembers[reducedChild].type == TYPE_DOUBLE_CHILD){
            if(nextChildMarker < 2){
                newCol->reducedMembers[reducedMember].childMarkerEdges[nextChildMarker] = markerOfParent(dec, findMember(dec,newCol->reducedMembers[reducedChild].member)); //TODO: check if find is necessary
                ++nextChildMarker;
            }
            newCol->reducedMembers[reducedMember].numTwoEnds++;
        }
    }
}

static void determineTypeParallel(
        MATRECGraphicColumnAddition *newCol, reduced_member_id reducedMemberId, int depth){
    const int numOneEnd = newCol->reducedMembers[reducedMemberId].numOneEnd;
    const int numTwoEnds = newCol->reducedMembers[reducedMemberId].numTwoEnds;
    assert(numOneEnd >= 0);
    assert(numTwoEnds >= 0);
    assert(numOneEnd + 2*numTwoEnds <= 2);
    MATRECColReducedMember * reducedMember = &newCol->reducedMembers[reducedMemberId];
    if( depth == 0){
        if(numTwoEnds + numOneEnd > 0 && pathEdgeIsValid(reducedMember->firstPathEdge)){
            reducedMember->type = TYPE_CYCLE_CHILD;
        }else{
            reducedMember->type = TYPE_ROOT;
        }
        return;
    }
    int numEnds = numOneEnd + 2*numTwoEnds;

    if(numEnds == 0 && pathEdgeIsValid(reducedMember->firstPathEdge)){
        reducedMember->type = TYPE_CYCLE_CHILD;
    }else if(numEnds == 1){
        reducedMember->type = TYPE_SINGLE_CHILD;
    }else if(numEnds == 2){
        if(pathEdgeIsValid(reducedMember->firstPathEdge)){
            newCol->remainsGraphic = false;
            reducedMember->type = TYPE_INVALID;
        }else{
            reducedMember->type = TYPE_DOUBLE_CHILD;
        }
    }else{
        //no child contains path edges, so we are a leaf of the reduced decomposition
        assert(pathEdgeIsValid(reducedMember->firstPathEdge));
        reducedMember->type = TYPE_CYCLE_CHILD; //TODO: is this not duplicate with first case? Should be able to turn into a switch case
    }
}
static void determineTypeSeries(MATRECGraphicDecomposition* dec, MATRECGraphicColumnAddition* newCol, reduced_member_id reducedMemberId,
                         int depth){
    const int numOneEnd = newCol->reducedMembers[reducedMemberId].numOneEnd;
    const int numTwoEnds = newCol->reducedMembers[reducedMemberId].numTwoEnds;
    assert(dec);
    assert(newCol);
    assert(numOneEnd >= 0);
    assert(numTwoEnds >= 0);
    assert(numOneEnd + 2*numTwoEnds <= 2);
    assert(getMemberType(dec, findMemberNoCompression(dec,newCol->reducedMembers[reducedMemberId].member)) == SPQR_MEMBERTYPE_SERIES);

    MATRECColReducedMember *reducedMember =&newCol->reducedMembers[reducedMemberId];
    spqr_member member = findMember(dec, reducedMember->member); //We could also pass this as function argument
    int countedPathEdges = 0;
    for(path_edge_id pathEdge = reducedMember->firstPathEdge; pathEdgeIsValid(pathEdge);
        pathEdge = newCol->pathEdges[pathEdge].nextMember){
        ++countedPathEdges;
    } //TODO: replace loop by count
    int numMemberEdges = getNumMemberEdges(dec,member);
    if(depth == 0){
        if(numTwoEnds != 0){
            if(countedPathEdges == numMemberEdges -1){
                reducedMember->type = TYPE_CYCLE_CHILD;
            }else{
                reducedMember->type = TYPE_INVALID;
                newCol->remainsGraphic = false;
            }

        }else{
            reducedMember->type = countedPathEdges == numMemberEdges-1 ? TYPE_CYCLE_CHILD : TYPE_ROOT;
        }
        return;
    }

    if(countedPathEdges == numMemberEdges-1){
        reducedMember->type = TYPE_CYCLE_CHILD;
    }else if(countedPathEdges + numTwoEnds == numMemberEdges-1){ //TODO: shouldn't this be numMemberEdges?
        assert(numTwoEnds == 1);
        reducedMember->type = TYPE_DOUBLE_CHILD;
    }else if(numTwoEnds == 1){
        newCol->remainsGraphic = false;
        reducedMember->type = TYPE_INVALID;
    }else{
        assert(numTwoEnds == 0);
        reducedMember->type = numOneEnd == 2 ? TYPE_DOUBLE_CHILD : TYPE_SINGLE_CHILD;
    }
}

static void determineTypeRigid(MATRECGraphicDecomposition* dec, MATRECGraphicColumnAddition* newCol, reduced_member_id reducedMemberId, int depth){
    //Rough idea; first, we find the
    const int numOneEnd = newCol->reducedMembers[reducedMemberId].numOneEnd;
    const int numTwoEnds = newCol->reducedMembers[reducedMemberId].numTwoEnds;
    assert(dec);
    assert(newCol);
    assert(numOneEnd >= 0);
    assert(numTwoEnds >= 0);
    assert(numOneEnd + 2*numTwoEnds <= 2);
    assert(getMemberType(dec, findMemberNoCompression(dec,newCol->reducedMembers[reducedMemberId].member)) == SPQR_MEMBERTYPE_RIGID);
    spqr_member member = findMember(dec, newCol->reducedMembers[reducedMemberId].member);

    spqr_node parentMarkerNodes[2] = {
            depth == 0 ? SPQR_INVALID_NODE : findEdgeHead(dec, markerToParent(dec, member)),
            depth == 0 ? SPQR_INVALID_NODE : findEdgeTail(dec, markerToParent(dec, member)),
    };

    spqr_edge * childMarkerEdges = newCol->reducedMembers[reducedMemberId].childMarkerEdges;
    spqr_node childMarkerNodes[4] = {
            childMarkerEdges[0] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeHead(dec, childMarkerEdges[0]),
            childMarkerEdges[0] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeTail(dec, childMarkerEdges[0]),
            childMarkerEdges[1] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeHead(dec, childMarkerEdges[1]),
            childMarkerEdges[1] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeTail(dec, childMarkerEdges[1]),
    };

    //First, find the end nodes of the path.
    //If there are too many (>4) or there is some node with degree > 2, we terminate
    spqr_node * pathEndNodes = newCol->reducedMembers[reducedMemberId].rigidEndNodes;
    for (int i = 0; i < 4; ++i) {
        pathEndNodes[i] = SPQR_INVALID_NODE;
    }
    int numPathEndNodes = 0;
    for (path_edge_id pathEdge = newCol->reducedMembers[reducedMemberId].firstPathEdge; pathEdgeIsValid(pathEdge);
         pathEdge = newCol->pathEdges[pathEdge].nextMember) {
        spqr_edge edge = newCol->pathEdges[pathEdge].edge;
        spqr_node nodes[2] = {findEdgeHead(dec, edge), findEdgeTail(dec, edge)};
        for (int i = 0; i < 2; ++i) {
            spqr_node node = nodes[i];
            assert(newCol->nodePathDegree[node] > 0);
            if(newCol->nodePathDegree[node] > 2){
                //Node of degree 3 or higher implies that the given edges cannot form a path.
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
            if(newCol->nodePathDegree[node] == 1){
                //If we have 5 or more end nodes, stop
                if(numPathEndNodes >= 4){
                    newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                    newCol->remainsGraphic = false;
                    return;
                }
                pathEndNodes[numPathEndNodes] = node;
                ++numPathEndNodes;
            }


        }
    }

    //Exchange end nodes in the end nodes so that if there are 4 edges, 0 is connected with 1 and 2 is connected with 3
    //Parent marker should come first for each path
    if(numPathEndNodes == 4){
        //We try to follow the path
        spqr_edge * nodeEdges;
        MATRECallocBlockArray(dec->env,&nodeEdges, (size_t) (2* largestNodeID(dec))); //TODO: move to struct

        //initialize for all relevant nodes
        for (path_edge_id pathEdge = newCol->reducedMembers[reducedMemberId].firstPathEdge; pathEdgeIsValid(pathEdge);
             pathEdge = newCol->pathEdges[pathEdge].nextMember){
            spqr_edge edge = newCol->pathEdges[pathEdge].edge;
            spqr_node nodes[2] = {findEdgeHead(dec, edge), findEdgeTail(dec, edge)};
            for (int i = 0; i < 2; ++i) {
                spqr_node node = nodes[i];
                nodeEdges[2*node] = SPQR_INVALID_EDGE;
                nodeEdges[2*node + 1] = SPQR_INVALID_EDGE;
            }
        }
        //store incident edges for every node
        for (path_edge_id pathEdge = newCol->reducedMembers[reducedMemberId].firstPathEdge; pathEdgeIsValid(pathEdge);
             pathEdge = newCol->pathEdges[pathEdge].nextMember){
            spqr_edge edge = newCol->pathEdges[pathEdge].edge;
            spqr_node nodes[2] = {findEdgeHead(dec, edge), findEdgeTail(dec, edge)};
            for (int i = 0; i < 2; ++i) {
                spqr_node node = nodes[i];
                spqr_node index = 2 * node;
                if(nodeEdges[index] != SPQR_INVALID_EDGE){
                    ++index;
                }
                nodeEdges[index] = edge;
            }
        }
        //Now follow the path starting from end node 0 to see where we end
        spqr_edge previousEdge = SPQR_INVALID_EDGE;
        spqr_node currentNode = pathEndNodes[0];
        while(true){
            spqr_edge edge = nodeEdges[2 * currentNode];
            if(edge == previousEdge){
                edge = nodeEdges[2*currentNode+1];
            }
            if(edge == SPQR_INVALID_EDGE){
                break;
            }
            previousEdge = edge;
            spqr_node node = findEdgeHead(dec, edge);
            currentNode = (node != currentNode) ? node : findEdgeTail(dec,edge);
        }
        MATRECfreeBlockArray(dec->env,&nodeEdges);

        if(currentNode == pathEndNodes[2]){
            pathEndNodes[2] = pathEndNodes[1];
            pathEndNodes[1] = currentNode;
        }else if(currentNode == pathEndNodes[3]){
            pathEndNodes[3] = pathEndNodes[1];
            pathEndNodes[1] = currentNode;
        }
        //make sure node 2 is the parent marker. We can assume that 2 and 3 now also form a nice path
        if(pathEndNodes[2] != parentMarkerNodes[0] && pathEndNodes[2] != parentMarkerNodes[1]){
            spqr_node temp = pathEndNodes[2];
            pathEndNodes[2] = pathEndNodes[3];
            pathEndNodes[3] = temp;
        }
    }
    //make sure node 0 is the parent marker node
    if(numPathEndNodes >= 2 && pathEndNodes[0] != parentMarkerNodes[0] && pathEndNodes[0] != parentMarkerNodes[1]){
        spqr_node temp = pathEndNodes[0];
        pathEndNodes[0] = pathEndNodes[1];
        pathEndNodes[1] = temp;
    }

    //Finally, we determine the type of the rigid node
    if(depth == 0){
        if(numPathEndNodes == 0){
            if(numOneEnd == 2 && (
                    childMarkerNodes[0] == childMarkerNodes[2] ||
                    childMarkerNodes[0] == childMarkerNodes[3] ||
                    childMarkerNodes[1] == childMarkerNodes[2] ||
                    childMarkerNodes[1] == childMarkerNodes[3])){
                newCol->reducedMembers[reducedMemberId].type = TYPE_ROOT;
                return;
            }else{
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
        }else if(numPathEndNodes == 2){
            if(numOneEnd == 1){
                bool pathAdjacentToChildMarker = false;
                bool childMarkerNodeAdjacent[2] = {false,false};
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        if(pathEndNodes[i] == childMarkerNodes[j]){
                            pathAdjacentToChildMarker = true;
                            childMarkerNodeAdjacent[j] = true;
                        }
                    }
                }

                if(childMarkerNodeAdjacent[0] && childMarkerNodeAdjacent[1]){
                    newCol->reducedMembers[reducedMemberId].type = TYPE_CYCLE_CHILD;
                    return;
                } else if(pathAdjacentToChildMarker){
                    newCol->reducedMembers[reducedMemberId].type = TYPE_ROOT;
                    return;
                }else{
                    newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                    newCol->remainsGraphic = false;
                    return;
                }
            }else if(numOneEnd == 2){
                bool childMarkerNodesMatched[2] = {false,false};
                bool endNodesMatched[2] ={false, false};
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        if(pathEndNodes[i] == childMarkerNodes[j]){
                            endNodesMatched[i] = true;
                            childMarkerNodesMatched[j/2] = true;
                        }
                    }
                }

                if(childMarkerNodesMatched[0] && childMarkerNodesMatched[1] && endNodesMatched[0] && endNodesMatched[1]){
                    newCol->reducedMembers[reducedMemberId].type = TYPE_ROOT;
                    return;
                }else{
                    newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                    newCol->remainsGraphic = false;
                    return;
                }
            }else if(numTwoEnds == 0){
                assert(numOneEnd == 0);
                newCol->reducedMembers[reducedMemberId].type = TYPE_ROOT;
                return;
            }else{
                assert(numOneEnd == 0);
                assert(numTwoEnds == 1);
                if((childMarkerNodes[0] == pathEndNodes[0] && childMarkerNodes[1] == pathEndNodes[1]) ||
                   (childMarkerNodes[0] == pathEndNodes[1] && childMarkerNodes[1] == pathEndNodes[0])
                        ){
                    newCol->reducedMembers[reducedMemberId].type = TYPE_CYCLE_CHILD;
                    return;
                }else{
                    newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                    newCol->remainsGraphic = false;
                    return;
                }
            }
        }else{
            assert(numPathEndNodes == 4);
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }
    }
    int parentMarkerDegrees[2] = {
            newCol->nodePathDegree[parentMarkerNodes[0]],
            newCol->nodePathDegree[parentMarkerNodes[1]]
    };
    //Non-root rigid member
    if(numPathEndNodes == 0){
        // We have no path edges, so there must be at least one child containing one/two path ends.
        assert(numOneEnd + numTwoEnds > 0);
        // We should not have a child marker edge parallel to the parent marker edge!
        assert(!(parentMarkerNodes[0] == childMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1])
               && !(parentMarkerNodes[0] == childMarkerNodes[1] && parentMarkerNodes[1] == childMarkerNodes[0]));
        if(numOneEnd == 0){
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }else if (numOneEnd == 1){
            if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
                || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1]){
                newCol->reducedMembers[reducedMemberId].type = TYPE_SINGLE_CHILD;
                return;
            }else{
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
        }else{
            assert(numOneEnd == 2);

            int childMarkerParentNode[2] = {-1,-1};
            bool isParallel = false;
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 2; ++j) {
                    if(childMarkerNodes[i] == parentMarkerNodes[j]){
                        if(childMarkerParentNode[i/2] >= 0){
                            isParallel = true;
                        }
                        childMarkerParentNode[i/2] = j;
                    }
                }
            }
            if(!isParallel && childMarkerParentNode[0] != -1 && childMarkerParentNode[1] != -1 &&
               childMarkerParentNode[0] != childMarkerParentNode[1]){
                newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
                return;
            }else{
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
        }
    }
    else if(numPathEndNodes == 2){
        if(numOneEnd == 1){
            //TODO: is the below check necessary?
            if(parentMarkerNodes[0] != pathEndNodes[0]){
                spqr_node tempMarker = parentMarkerNodes[0];
                parentMarkerNodes[0] = parentMarkerNodes[1];
                parentMarkerNodes[1] = tempMarker;

                int tempDegree = parentMarkerDegrees[0];
                parentMarkerDegrees[0] = parentMarkerDegrees[1];
                parentMarkerDegrees[1] = tempDegree;
            }
            if(parentMarkerNodes[0] != pathEndNodes[0]){
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
            if(parentMarkerNodes[1] == pathEndNodes[1]){
                // Path closes a cycle with parent marker edge.
                if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
                    || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1])
                {
                    newCol->reducedMembers[reducedMemberId].type = TYPE_SINGLE_CHILD;
                    return;
                }
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }else{
                if(childMarkerNodes[0] == pathEndNodes[1] || childMarkerNodes[1] == pathEndNodes[1]){
                    newCol->reducedMembers[reducedMemberId].type = TYPE_SINGLE_CHILD;
                    return;
                }else if(childMarkerNodes[0] == parentMarkerNodes[1] || childMarkerNodes[1] == parentMarkerNodes[1]){
                    newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
                    return;
                }
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
        }
        else if(numOneEnd == 2){
            spqr_node otherParentNode;
            if(pathEndNodes[0] == parentMarkerNodes[0]){
                otherParentNode = parentMarkerNodes[1];
            }else if(pathEndNodes[0] == parentMarkerNodes[1]){
                otherParentNode = parentMarkerNodes[0];
            }else{
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
            if(pathEndNodes[1] == otherParentNode){
                newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
                newCol->remainsGraphic = false;
                return;
            }
            bool childMatched[2] = {false,false};
            bool pathEndMatched = false;
            bool otherParentMatched = false;
            for (int i = 0; i < 4; ++i) {
                if(childMarkerNodes[i] == pathEndNodes[1]){
                    childMatched[i/2] = true;
                    pathEndMatched = true;
                }
                if(childMarkerNodes[i] == otherParentNode){
                    childMatched[i/2] = true;
                    otherParentMatched = true;
                }
            }
            if(childMatched[0] && childMatched[1] && pathEndMatched  && otherParentMatched){
                newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
                return;
            }
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;

        }
        else if(numTwoEnds == 0){
            if ((parentMarkerDegrees[0] % 2 == 0 && parentMarkerDegrees[1] == 1) ||
                (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] % 2 == 0))
            {
                newCol->reducedMembers[reducedMemberId].type = TYPE_SINGLE_CHILD;
                return;
            }
            else if (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] == 1)
            {
                newCol->reducedMembers[reducedMemberId].type = TYPE_CYCLE_CHILD;
                return;
            }
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;

        }
        else{
            assert(numTwoEnds == 1);
            if ((pathEndNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[0]
                 && childMarkerNodes[1] == pathEndNodes[1])
                || (pathEndNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1]
                    && childMarkerNodes[0] == pathEndNodes[1])
                || (pathEndNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[0]
                    && childMarkerNodes[1] == pathEndNodes[1])
                || (pathEndNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[1]
                    && childMarkerNodes[0] == pathEndNodes[1]))
            {
                newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
                return;
            }
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }
    }
    else {
        assert(numPathEndNodes == 4);
        if(pathEndNodes[0] != parentMarkerNodes[0] && pathEndNodes[0] != parentMarkerNodes[1]){
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }
        if(pathEndNodes[2] != parentMarkerNodes[0] && pathEndNodes[2] != parentMarkerNodes[1]){
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }
        if(numOneEnd == 1){
            if((pathEndNodes[1] == childMarkerNodes[0] || pathEndNodes[1] == childMarkerNodes[1]) ||
               (pathEndNodes[3] == childMarkerNodes[0] || pathEndNodes[3] == childMarkerNodes[1])){
                newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
                return;
            }
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }
        else if(numOneEnd == 2){
            bool pathConnected[2] = {false,false};
            bool childConnected[2] = {false,false};
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 4; ++j) {
                    if(pathEndNodes[1+2*i] == childMarkerNodes[j]){
                        pathConnected[i] = true;
                        childConnected[j/2] = true;
                    }
                }
            }
            if(pathConnected[0] && pathConnected[1] && childConnected[0] && childConnected[1]){
                newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
                return;
            }
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }
        else if(numTwoEnds == 0){
            newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
            return;
        }
        else{
            if((pathEndNodes[1] == childMarkerNodes[0] && pathEndNodes[3] == childMarkerNodes[1]) ||
               (pathEndNodes[1] == childMarkerNodes[1] && pathEndNodes[3] == childMarkerNodes[0])){
                newCol->reducedMembers[reducedMemberId].type = TYPE_DOUBLE_CHILD;
                return;
            }
            newCol->reducedMembers[reducedMemberId].type = TYPE_INVALID;
            newCol->remainsGraphic = false;
            return;
        }
    }
}
static void determineTypes(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol, MATRECColReducedComponent * component,
                    reduced_member_id reducedMember,
                    int depth ){
    assert(dec);
    assert(newCol);

    for (children_idx idx = newCol->reducedMembers[reducedMember].firstChild;
         idx < newCol->reducedMembers[reducedMember].firstChild
               + newCol->reducedMembers[reducedMember].numChildren;
         ++idx) {
        reduced_member_id reducedChild = newCol->childrenStorage[idx];
        assert(reducedMemberIsValid(reducedChild));
        determineTypes(dec,newCol,component,reducedChild,depth + 1);
        if(!newCol->remainsGraphic){
            return;
        }
    }

    countChildrenTypes(dec,newCol,reducedMember);
    if(2*newCol->reducedMembers[reducedMember].numTwoEnds + newCol->reducedMembers[reducedMember].numOneEnd > 2){
        newCol->remainsGraphic = false;
        return;
    }
    //Determine type of this
    bool isRoot = reducedMember == component->root;
    spqr_member member = findMember(dec, newCol->reducedMembers[reducedMember].member); //TODO: find necessary?
    SPQRMemberType type = getMemberType(dec, member);
    if(type == SPQR_MEMBERTYPE_PARALLEL){
        determineTypeParallel(newCol,reducedMember,depth);
    }else if (type == SPQR_MEMBERTYPE_SERIES){
        determineTypeSeries(dec,newCol,reducedMember,depth);
    }else if (type == SPQR_MEMBERTYPE_LOOP){
        newCol->reducedMembers[reducedMember].type = TYPE_ROOT;
    } else{
        assert(type == SPQR_MEMBERTYPE_RIGID);
        determineTypeRigid(dec,newCol,reducedMember,depth);
    }

    //Add a marked edge to the path edge of the parent of this
    if(newCol->remainsGraphic && !isRoot && newCol->reducedMembers[reducedMember].type == TYPE_CYCLE_CHILD){
        spqr_member parentMember = findMemberParent(dec, newCol->reducedMembers[reducedMember].member);
        reduced_member_id reducedParent = newCol->memberInformation[parentMember].reducedMember;
        spqr_edge marker = markerOfParent(dec, member);

        createPathEdge(dec,newCol,marker,reducedParent);
    }
}
static void propagateRootCycle(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition * newCol,MATRECColReducedComponent * component){
    reduced_member_id root = component->root;
    reduced_member_id uniqueNonPropagatedChild = INVALID_REDUCED_MEMBER;

    for (children_idx idx = newCol->reducedMembers[root].firstChild;
         idx < newCol->reducedMembers[root].firstChild + newCol->reducedMembers[root].numChildren; ++idx) {
        reduced_member_id reducedChild = newCol->childrenStorage[idx];
        if(newCol->reducedMembers[reducedChild].type != TYPE_CYCLE_CHILD){
            if(reducedMemberIsValid(uniqueNonPropagatedChild)){
                uniqueNonPropagatedChild = INVALID_REDUCED_MEMBER;
                break;
            }
            uniqueNonPropagatedChild = reducedChild;
        }
    }

    while(reducedMemberIsValid(uniqueNonPropagatedChild) && newCol->reducedMembers[root].type == TYPE_CYCLE_CHILD){
        spqr_edge edge = markerToParent(dec, findMember(dec,newCol->reducedMembers[uniqueNonPropagatedChild].member));
        createPathEdge(dec,newCol,edge,uniqueNonPropagatedChild);

        component->root = uniqueNonPropagatedChild;
        ++component->rootDepth;
        root = uniqueNonPropagatedChild;
        spqr_member member = findMember(dec, newCol->reducedMembers[root].member); //TODO: find necessary?
        SPQRMemberType type = getMemberType(dec, member);
        if(type == SPQR_MEMBERTYPE_PARALLEL){
            determineTypeParallel(newCol,root,0);
        }else if (type == SPQR_MEMBERTYPE_SERIES){
            determineTypeSeries(dec,newCol,root,0);
        }else{
            assert(type == SPQR_MEMBERTYPE_RIGID);
            determineTypeRigid(dec,newCol,root,0);
        }
        if(!newCol->remainsGraphic){
            return;
        }

        uniqueNonPropagatedChild = INVALID_REDUCED_MEMBER;
        for (children_idx idx = newCol->reducedMembers[root].firstChild;
             idx < newCol->reducedMembers[root].firstChild + newCol->reducedMembers[root].numChildren; ++idx) {
            reduced_member_id reducedChild = newCol->childrenStorage[idx];
            if(newCol->reducedMembers[reducedChild].type != TYPE_CYCLE_CHILD){
                if(reducedMemberIsValid(uniqueNonPropagatedChild)){
                    uniqueNonPropagatedChild = INVALID_REDUCED_MEMBER;
                    break;
                }
                uniqueNonPropagatedChild = reducedChild;
            }
        }
        if(type == SPQR_MEMBERTYPE_PARALLEL && reducedMemberIsValid(uniqueNonPropagatedChild)){
            newCol->reducedMembers[root].type = TYPE_CYCLE_CHILD;
        }else if (type == SPQR_MEMBERTYPE_RIGID && reducedMemberIsValid(uniqueNonPropagatedChild)){
            spqr_edge rigidMarker = markerOfParent(dec, findMember(dec,newCol->reducedMembers[uniqueNonPropagatedChild].member));
            assert(SPQRmemberIsValid(findEdgeMemberNoCompression(dec,rigidMarker)));
            int numEndNodes = 0;
            for (int i = 0; i < 4; ++i) {
                if(SPQRnodeIsInvalid(newCol->reducedMembers[root].rigidEndNodes[i])){
                   break;
                }
                ++numEndNodes;
            }
            if(numEndNodes == 2){
                spqr_node * pathEndNodes = newCol->reducedMembers[root].rigidEndNodes;
                spqr_node markerNodes[2] = {findEdgeTail(dec,rigidMarker),findEdgeHead(dec,rigidMarker)};
                if((markerNodes[0] == pathEndNodes[0] && markerNodes[1] == pathEndNodes[1]) ||
                (markerNodes[0] == pathEndNodes[1] && markerNodes[1] == pathEndNodes[0])){
                    newCol->reducedMembers[root].type = TYPE_CYCLE_CHILD;
                }
            }
        }
    }
}
static void determineComponentTypes(MATRECGraphicDecomposition * dec, MATRECGraphicColumnAddition * newCol, MATRECColReducedComponent * component){
    assert(dec);
    assert(newCol);
    assert(component);

    determineTypes(dec,newCol,component,component->root,0);
    if(newCol->remainsGraphic){
        propagateRootCycle(dec,newCol,component);
    }
}

MATREC_ERROR
MATRECGraphicColumnAdditionCheck(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol, MATREC_col column, const MATREC_row *rows, size_t numRows) {
    assert(dec);
    assert(newCol);
    assert(numRows == 0 || rows);

    newCol->remainsGraphic = true;
    cleanupPreviousIteration(dec, newCol);
    //assert that previous iteration was cleaned up

    //Store call data
    MATREC_CALL(newColUpdateColInformation(dec, newCol, column, rows, numRows));

    //compute reduced decomposition
    MATREC_CALL(constructReducedDecomposition(dec, newCol));
    //initialize path edges in reduced decomposition
    MATREC_CALL(createPathEdges(dec,newCol));
    //determine types
    for (int i = 0; i < newCol->numReducedComponents; ++i) {
        determineComponentTypes(dec,newCol,&newCol->reducedComponents[i]);
    }
    //clean up memberInformation
    cleanUpMemberInformation(newCol);

    return MATREC_OKAY;
}

///Contains the data which tells us where to store the new column after the graph has been modified
///In case member is a parallel or series node, the respective new column and rows are placed in parallel (or series) with it
///Otherwise, the rigid member has a free spot between firstNode and secondNode
typedef struct {
    spqr_member member;
    spqr_node terminalNode[2];
    int numTerminals;
} NewColInformation;

static NewColInformation emptyNewColInformation(void){
    NewColInformation information;
    information.member = SPQR_INVALID_MEMBER;
    information.terminalNode[0] = SPQR_INVALID_NODE;
    information.terminalNode[1] = SPQR_INVALID_NODE;
    information.numTerminals = 0;
    return information;
}
static void addTerminal(NewColInformation * info, spqr_node node){
    assert(info->numTerminals < 2);
    info->terminalNode[info->numTerminals] = node;
    info->numTerminals += 1;
}
static void setTerminalMember(NewColInformation * info, spqr_member member){
    info->member = member;
}

static MATREC_ERROR mergeMemberIntoParent(MATRECGraphicDecomposition *dec, spqr_member member,
                                        bool headToHead // controls the orientation of the merge, e.g. merge heads of the two edges into node if headToHead is true
){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));
    spqr_member parentMember = findMemberParent(dec, member);

    assert(SPQRmemberIsValid(parentMember));
    assert(memberIsRepresentative(dec,parentMember));

    spqr_edge parentToChild = markerOfParent(dec, member);
    removeEdgeFromMemberEdgeList(dec,parentToChild,parentMember);
    spqr_edge childToParent = markerToParent(dec, member);
    removeEdgeFromMemberEdgeList(dec,childToParent,member);

    spqr_node parentEdgeNodes[2] = {findEdgeTail(dec, parentToChild), findEdgeHead(dec, parentToChild)};
    spqr_node childEdgeNodes[2] = {findEdgeTail(dec, childToParent), findEdgeHead(dec, childToParent)};

    clearEdgeHeadAndTail(dec,parentToChild);
    clearEdgeHeadAndTail(dec,childToParent);

    spqr_node first = childEdgeNodes[headToHead ? 0 : 1];
    spqr_node second = childEdgeNodes[headToHead ? 1 : 0];
    {
        spqr_node newNode = mergeNodes(dec, parentEdgeNodes[0], first);
        spqr_node toRemoveFrom = newNode == first ? parentEdgeNodes[0] : first;
        mergeNodeEdgeList(dec,newNode,toRemoveFrom);
    }
    {
        spqr_node newNode = mergeNodes(dec, parentEdgeNodes[1], second);
        spqr_node toRemoveFrom = newNode == second ? parentEdgeNodes[1] : second;
        mergeNodeEdgeList(dec,newNode,toRemoveFrom);
    }


    spqr_member newMember = mergeMembers(dec, member, parentMember);
    spqr_member toRemoveFrom = newMember == member ? parentMember : member;
    mergeMemberEdgeList(dec,newMember,toRemoveFrom);
    if(toRemoveFrom == parentMember){
        updateMemberParentInformation(dec,newMember,toRemoveFrom);
    }
    updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);

    return MATREC_OKAY;
}



static MATREC_ERROR createParallelNodes(MATRECGraphicDecomposition * dec, spqr_member member){
    spqr_node head,tail;
    MATREC_CALL(createNode(dec,&head));
    MATREC_CALL(createNode(dec,&tail));

    spqr_edge firstEdge = getFirstMemberEdge(dec, member);
    spqr_edge edge = firstEdge;
    do{
        setEdgeHeadAndTail(dec,edge,head,tail);
        edge = getNextMemberEdge(dec,edge);
    }while(edge != firstEdge);
    return MATREC_OKAY;
}

static MATREC_ERROR splitParallel(MATRECGraphicDecomposition *dec, spqr_member parallel,
                                spqr_edge edge1, spqr_edge edge2,
                                spqr_member * childParallel){
    assert(dec);
    assert(SPQRedgeIsValid(edge1));
    assert(SPQRedgeIsValid(edge2));
    assert(SPQRmemberIsValid(parallel));

    bool childContainsTree = edgeIsTree(dec,edge1) || edgeIsTree(dec,edge2);
    spqr_edge toParent = markerToParent(dec, parallel);
    bool parentMoved = toParent == edge1 || toParent == edge2;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, childParallel));

    moveEdgeToNewMember(dec,edge1,parallel,*childParallel);
    moveEdgeToNewMember(dec,edge2,parallel,*childParallel);

    if(parentMoved){
        MATREC_CALL(createMarkerPair(dec,*childParallel,parallel,!childContainsTree));
    }else{
        MATREC_CALL(createMarkerPair(dec,parallel,*childParallel,childContainsTree));
    }
    return MATREC_OKAY;
}
static MATREC_ERROR transformParallel(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol, reduced_member_id reducedMemberId,
                                    NewColInformation * newColInfo, int depth){
    assert(dec);
    assert(newCol);
    MATRECColReducedMember * reducedMember = &newCol->reducedMembers[reducedMemberId];
    if(depth != 0){
        assert(reducedMember->numOneEnd == 1);
        //TODO: split parallel
        if(getNumMemberEdges(dec,reducedMember->member) > 3){
            spqr_member child;
            MATREC_CALL(splitParallel(dec,reducedMember->member,reducedMember->childMarkerEdges[0],
                                    markerToParent(dec,reducedMember->member),&child));
            assert(memberIsRepresentative(dec,child));
            reducedMember->member = child;
        }
        //TODO: can cleanup parallel merging to not immediately identify the nodes by assigning instead
        MATREC_CALL(createParallelNodes(dec,reducedMember->member));
        MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec,reducedMember->childMarkerEdges[0]),
                                        pathEdgeIsInvalid(reducedMember->firstPathEdge)));
        return MATREC_OKAY;
    }
    if(reducedMember->numOneEnd == 0 && reducedMember->numTwoEnds == 0){
        assert(pathEdgeIsValid(reducedMember->firstPathEdge) && reducedMember->numPathEdges == 1);
        //The new edge can be placed in parallel; just add it to this member
        setTerminalMember(newColInfo,reducedMember->member);
        return MATREC_OKAY;
    }
    assert(reducedMember->numOneEnd == 2);

    //split off if the parallel contains more than three edges
    if(getNumMemberEdges(dec,reducedMember->member) > 3){
        spqr_member child;
        MATREC_CALL(splitParallel(dec,reducedMember->member,reducedMember->childMarkerEdges[0],
                                reducedMember->childMarkerEdges[1],&child));
        assert(memberIsRepresentative(dec,child));
        reducedMember->member = child;
    }
    assert(getNumMemberEdges(dec,reducedMember->member) == 3);

    //TODO: again, this can probably be done more efficiently
    MATREC_CALL(createParallelNodes(dec,reducedMember->member));
    MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec,reducedMember->childMarkerEdges[0]),true));

    //If this doesn't work, try;
    bool headToHead = reducedMember->numPathEdges == 0;
//    bool headToHead = reducedMember->numPathEdges == 0 && reducedMember->type != TYPE_DOUBLE_CHILD && reducedMember->type != TYPE_ROOT;
    MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec,reducedMember->childMarkerEdges[1]),headToHead));
    setTerminalMember(newColInfo, findMember(dec,reducedMember->member));
    return MATREC_OKAY;
}

static MATREC_ERROR splitSeries(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition * newCol,
                              MATRECColReducedMember * reducedMember,
                              spqr_member member, spqr_member * loopMember){
    assert(dec);
    assert(reducedMember);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));

    bool createPathSeries = reducedMember->numPathEdges > 1;
    bool convertOriginal = reducedMember->numPathEdges == getNumMemberEdges(dec,member) -1;
    if(!createPathSeries && convertOriginal){
        //only one path edge; we are in a loop; no need to change anything
        assert(getNumMemberEdges(dec,member) == 2);
        assert(reducedMember->numPathEdges == 1);
        *loopMember = member;
        changeLoopToParallel(dec,member);
        return MATREC_OKAY;
    }

    spqr_member pathMember;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

    path_edge_id pathEdgeId = reducedMember->firstPathEdge;
    bool parentMoved = false;
    while(pathEdgeIsValid(pathEdgeId)){
        spqr_edge pathEdge = newCol->pathEdges[pathEdgeId].edge;
        pathEdgeId = newCol->pathEdges[pathEdgeId].nextMember;
        if(pathEdge == markerToParent(dec,member)){
            parentMoved = true;
        }
        moveEdgeToNewMember(dec,pathEdge,member,pathMember);
    }
    if(convertOriginal == createPathSeries){
        if(parentMoved){
            MATREC_CALL(createMarkerPair(dec,pathMember,member,false));
        }else{
            MATREC_CALL(createMarkerPair(dec,member,pathMember,true));
        }
        *loopMember = convertOriginal ? member : pathMember;
        changeLoopToParallel(dec,*loopMember);
        return MATREC_OKAY;
    }
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, loopMember));
    if(parentMoved){
        MATREC_CALL(createMarkerPair(dec,pathMember,*loopMember,false));
        MATREC_CALL(createMarkerPair(dec,*loopMember,member,false));
    }else{
        MATREC_CALL(createMarkerPair(dec,member,*loopMember,true));
        MATREC_CALL(createMarkerPair(dec,*loopMember,pathMember,true));
    }

    return MATREC_OKAY;
}


static MATREC_ERROR splitSeriesMerging(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition * newCol,
                                     MATRECColReducedMember * reducedMember,
                                     spqr_member member,
                                     spqr_edge * pathRepresentative,
                                     spqr_edge * nonPathRepresentative,
                                     spqr_edge exceptionEdge1,
                                     spqr_edge exceptionEdge2){
    assert(dec);
    assert(reducedMember);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));

    int numExceptionEdges = (exceptionEdge1 == SPQR_INVALID_EDGE ? 0 : 1) + (exceptionEdge2 == SPQR_INVALID_EDGE ? 0 : 1);
    int numNonPathEdges = getNumMemberEdges(dec,member) - reducedMember->numPathEdges - numExceptionEdges;
    bool createPathSeries = reducedMember->numPathEdges > 1;
    //If this holds, there are 2 or more non-parent marker non-path edges
    bool createNonPathSeries = numNonPathEdges > 1;
    assert(exceptionEdge1 == SPQR_INVALID_EDGE || !newCol->edgeInPath[exceptionEdge1]);
    assert(exceptionEdge2 == SPQR_INVALID_EDGE || !newCol->edgeInPath[exceptionEdge2]);

    if(createPathSeries){
        spqr_member pathMember;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

        path_edge_id pathEdgeId = reducedMember->firstPathEdge;
        bool parentMoved = false;
        while(pathEdgeIsValid(pathEdgeId)){
            spqr_edge pathEdge = newCol->pathEdges[pathEdgeId].edge;
            pathEdgeId = newCol->pathEdges[pathEdgeId].nextMember;
            assert(pathEdge != exceptionEdge1 && pathEdge != exceptionEdge2);
            parentMoved = parentMoved || markerToParent(dec,member) == pathEdge;
            moveEdgeToNewMember(dec,pathEdge,member,pathMember);
        }
        assert(getNumMemberEdges(dec,pathMember) >= 2);

        spqr_edge ignored;
        if(parentMoved){
            MATREC_CALL(createMarkerPairWithReferences(dec,pathMember,member,false,&ignored,pathRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,pathMember,true,pathRepresentative,&ignored));
        }
    }else{
        if(pathEdgeIsValid(reducedMember->firstPathEdge)){
            *pathRepresentative = newCol->pathEdges[reducedMember->firstPathEdge].edge;
        }else{
            *pathRepresentative = SPQR_INVALID_EDGE;
        }
    }

    if(createNonPathSeries){
        spqr_member nonPathMember;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &nonPathMember));

        spqr_edge edge = getFirstMemberEdge(dec, member);
        bool parentMoved = false;
        bool canStop = false; //hack when the first edge is moved in the below loop to prevent that we immediately terminate
        do{
            spqr_edge nextEdge = getNextMemberEdge(dec, edge);
            if(edge != *pathRepresentative && edge != exceptionEdge1 && edge != exceptionEdge2){
                parentMoved = parentMoved || markerToParent(dec,member) == edge;
                moveEdgeToNewMember(dec,edge,member,nonPathMember);
            }else{
                canStop = true;
            }
            edge = nextEdge;
            if(canStop && edge == getFirstMemberEdge(dec,member)){
                break;
            }
        }while(true);
        assert(getNumMemberEdges(dec,nonPathMember) >= 2);
        bool representativeIsTree = !edgeIsTree(dec,exceptionEdge1);
        if(SPQRedgeIsValid(exceptionEdge2)){
            representativeIsTree = representativeIsTree || !edgeIsTree(dec,exceptionEdge2);
        }
        spqr_edge ignored;
        if(parentMoved){
            MATREC_CALL(createMarkerPairWithReferences(dec,nonPathMember,member,!representativeIsTree,&ignored,nonPathRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,nonPathMember,representativeIsTree,nonPathRepresentative,&ignored));
        }
    }else{
        *nonPathRepresentative = SPQR_INVALID_EDGE;
        if(numNonPathEdges != 0) {
            spqr_edge firstEdge = getFirstMemberEdge(dec, member);
            spqr_edge edge = firstEdge;
            do {
                if (edge != *pathRepresentative && edge != exceptionEdge1 && edge != exceptionEdge2) {
                    *nonPathRepresentative = edge;
                    break;
                }
                edge = getNextMemberEdge(dec, edge);
            } while (edge != firstEdge);
            assert(*nonPathRepresentative != SPQR_INVALID_EDGE);
        }
    }

    return MATREC_OKAY;
}

static MATREC_ERROR transformSeries(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol, reduced_member_id reducedMemberId,
                                  NewColInformation * newColInfo, int depth){
    assert(dec);
    assert(newCol);
    assert(newColInfo);
    MATRECColReducedMember * reducedMember = &newCol->reducedMembers[reducedMemberId];
    if(depth == 0){
        // If we have a child containing both ends then we should have moved the reduced root there.
        //This means that in this case, we know this is the only reduced member
        assert(reducedMember->numTwoEnds == 0);
        if(reducedMember->numPathEdges == getNumMemberEdges(dec,reducedMember->member)-1){
            spqr_member adjacentMember = SPQR_INVALID_MEMBER;
            {
                spqr_edge firstEdge = getFirstMemberEdge(dec, reducedMember->member);
                spqr_edge edge = firstEdge;
                do{
                    if(!newCol->edgeInPath[edge]){
                        if(edge == markerToParent(dec,reducedMember->member)){
                            adjacentMember = findMemberParent(dec,reducedMember->member);
                        }else if(edgeIsMarker(dec,edge)){
                            adjacentMember = findEdgeChildMember(dec,edge);
                        }

                        break; //There is only a singular such edge
                    }
                    edge = getNextMemberEdge(dec,edge);
                }while(edge != firstEdge);
            }

            if(SPQRmemberIsValid(adjacentMember)){
                SPQRMemberType adjacentType = getMemberType(dec, adjacentMember);
                if(adjacentType == SPQR_MEMBERTYPE_PARALLEL){
                    setTerminalMember(newColInfo,adjacentMember);
                    return MATREC_OKAY;
                }
            }
        }

        if(reducedMember->numOneEnd == 0){
            //Isolated single cycle
            spqr_member loopMember;
            MATREC_CALL(splitSeries(dec,newCol,reducedMember,reducedMember->member,&loopMember));
            setTerminalMember(newColInfo,loopMember);
            return MATREC_OKAY;
        }
        if(reducedMember->numOneEnd == 1){
            //Cycle is root of multiple components
            //Split off all path edges and all non-path edges, except the unique child marker
            spqr_edge childMarker = reducedMember->childMarkerEdges[0];
            assert(childMarker != SPQR_INVALID_EDGE);

            spqr_edge pathRepresentative,nonPathRepresentative;
            MATREC_CALL(splitSeriesMerging(dec, newCol, reducedMember, reducedMember->member, &pathRepresentative, &nonPathRepresentative,
                                         childMarker, SPQR_INVALID_EDGE));
            assert(getNumMemberEdges(dec,reducedMember->member) == 3);

            spqr_node a,b,c;
            MATREC_CALL(createNode(dec,&a));
            MATREC_CALL(createNode(dec,&b));
            MATREC_CALL(createNode(dec,&c));

            assert(pathRepresentative != childMarker &&
                   nonPathRepresentative != childMarker &&
                   pathRepresentative != nonPathRepresentative);
            assert(edgeIsTree(dec,pathRepresentative) );

            setEdgeHeadAndTail(dec,childMarker,a,b);
            setEdgeHeadAndTail(dec,pathRepresentative,b,c);
            setEdgeHeadAndTail(dec,nonPathRepresentative,c,a);

            addTerminal(newColInfo,c);

            //parent path always ends in the tail node, so we identify head to head
            MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec,childMarker),true));

            spqr_member member = findMember(dec, reducedMember->member);
            setTerminalMember(newColInfo,member);

            return MATREC_OKAY;
        }
        assert(reducedMember->numOneEnd == 2);
        // If there is more than 1 (non)-path edge, we split off by moving them to a new series member and creating a parallel
        spqr_edge pathEdge = SPQR_INVALID_EDGE;
        spqr_edge nonPathEdge = SPQR_INVALID_EDGE;

        MATREC_CALL(splitSeriesMerging(dec, newCol, reducedMember, reducedMember->member, &pathEdge, &nonPathEdge,
                                     reducedMember->childMarkerEdges[0],
                                     reducedMember->childMarkerEdges[1]));

        //Check that all edges are unique
        assert(reducedMember->childMarkerEdges[0] != reducedMember->childMarkerEdges[1]);
        assert(pathEdge != reducedMember->childMarkerEdges[0] && pathEdge != reducedMember->childMarkerEdges[1] &&
               nonPathEdge != reducedMember->childMarkerEdges[0] && nonPathEdge != reducedMember->childMarkerEdges[1]);
        assert(pathEdge != SPQR_INVALID_EDGE || nonPathEdge != SPQR_INVALID_EDGE);
        assert(getNumMemberEdges(dec, reducedMember->member) == 3 ||
               getNumMemberEdges(dec, reducedMember->member) == 4);
        //After splitting there is the following possibilities for nodes a-d:
        //(a)-child-(b)-path-(c)-child-(d)-nonpath-(a)
        //(a)-child-(b)-path-(c)-child-(d==a)
        //(a)-child-(b)=(c)-child-(d)-nonpath-(a)

        spqr_node a = SPQR_INVALID_NODE;
        spqr_node b = SPQR_INVALID_NODE;
        spqr_node c = SPQR_INVALID_NODE;
        spqr_node d = SPQR_INVALID_NODE;
        MATREC_CALL(createNode(dec, &a));
        MATREC_CALL(createNode(dec, &b));
        if (SPQRedgeIsValid(pathEdge)) {
            MATREC_CALL(createNode(dec, &c));
        } else {
            c = b;
        }
        if (SPQRedgeIsValid(nonPathEdge)) {
            MATREC_CALL(createNode(dec, &d));
        } else {
            d = a;
        }

        setEdgeHeadAndTail(dec, reducedMember->childMarkerEdges[0], a, b);
        setEdgeHeadAndTail(dec, reducedMember->childMarkerEdges[1], d, c);
        if (SPQRedgeIsValid(pathEdge)) {
            setEdgeHeadAndTail(dec, pathEdge, b, c);
        }
        if (SPQRedgeIsValid(nonPathEdge)) {
            setEdgeHeadAndTail(dec, nonPathEdge, d, a);
        }
        MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec, reducedMember->childMarkerEdges[0]), true));
        MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec, reducedMember->childMarkerEdges[1]), true));
        setTerminalMember(newColInfo, findMember(dec,reducedMember->member));

        return MATREC_OKAY;
    }
    assert(reducedMember->type == TYPE_SINGLE_CHILD);
    spqr_edge parentMarker = markerToParent(dec, reducedMember->member);
    assert(SPQRedgeIsValid(parentMarker));

    if (reducedMember->numOneEnd == 0) {
        assert(reducedMember->numTwoEnds == 0);
        assert(pathEdgeIsValid(reducedMember->firstPathEdge));
        //propagate path edge

        //Split off all path edges and all non-path edges EXCEPT the parent marker if there are more than 2
        spqr_edge pathRepresentative, nonPathRepresentative;
        MATREC_CALL(splitSeriesMerging(dec, newCol, reducedMember, reducedMember->member, &pathRepresentative,
                                     &nonPathRepresentative,
                                     parentMarker, SPQR_INVALID_EDGE));

        assert(getNumMemberEdges(dec, reducedMember->member) == 3);

        //add nodes edges and a single terminal

        spqr_node a, b, c;
        MATREC_CALL(createNode(dec, &a));
        MATREC_CALL(createNode(dec, &b));
        MATREC_CALL(createNode(dec, &c));

        assert(pathRepresentative != parentMarker && nonPathRepresentative != parentMarker &&
               pathRepresentative != nonPathRepresentative);

        setEdgeHeadAndTail(dec, parentMarker, a, b);
        setEdgeHeadAndTail(dec, pathRepresentative, b, c);
        setEdgeHeadAndTail(dec, nonPathRepresentative, c, a);

        addTerminal(newColInfo, c);
        return MATREC_OKAY;
    }

    //Path passes through this marker and the parent marker
    assert(reducedMember->numOneEnd == 1);

    spqr_edge childMarker = reducedMember->childMarkerEdges[0];
    assert(SPQRedgeIsValid(childMarker));

    // If there is more than 1 (non)-path edge, we split off by moving them to a new series member and creating a parallel
    spqr_edge pathEdge = SPQR_INVALID_EDGE;
    spqr_edge nonPathEdge = SPQR_INVALID_EDGE;

    MATREC_CALL(splitSeriesMerging(dec, newCol, reducedMember, reducedMember->member, &pathEdge, &nonPathEdge,
                                 childMarker,parentMarker));

    //Check that all edges are unique
    assert(childMarker != parentMarker);
    assert(pathEdge != childMarker && pathEdge != parentMarker &&
           nonPathEdge != childMarker && nonPathEdge != parentMarker);
    assert(pathEdge != SPQR_INVALID_EDGE || nonPathEdge != SPQR_INVALID_EDGE);
    assert(getNumMemberEdges(dec, reducedMember->member) == 3 ||
           getNumMemberEdges(dec, reducedMember->member) == 4);
    //After splitting there is the following possibilities for nodes a-d:
    //(a)-parent-(b)-path-(c)-child-(d)-nonpath-(a)
    //(a)-parent-(b)-path-(c)-child-(d==a)
    //(a)-parent-(b)=(c)-child-(d)-nonpath-(a)

    spqr_node a = SPQR_INVALID_NODE;
    spqr_node b = SPQR_INVALID_NODE;
    spqr_node c = SPQR_INVALID_NODE;
    spqr_node d = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec, &a));
    MATREC_CALL(createNode(dec, &b));
    if (SPQRedgeIsValid(pathEdge)) {
        MATREC_CALL(createNode(dec, &c));
    } else {
        c = b;
    }
    if (SPQRedgeIsValid(nonPathEdge)) {
        MATREC_CALL(createNode(dec, &d));
    } else {
        d = a;
    }

    setEdgeHeadAndTail(dec, parentMarker, a, b);
    setEdgeHeadAndTail(dec, childMarker, d, c);
    if (SPQRedgeIsValid(pathEdge)) {
        setEdgeHeadAndTail(dec, pathEdge, b, c);
    }
    if (SPQRedgeIsValid(nonPathEdge)) {
        setEdgeHeadAndTail(dec, nonPathEdge, d, a);
    }
    MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec, childMarker), true));

    return MATREC_OKAY;
}
static MATREC_ERROR transformRigid(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition * newCol, reduced_member_id id,
                                 NewColInformation * newColInfo, int depth){
    assert(dec);
    assert(newCol);
    assert(reducedMemberIsValid(id));
    assert(depth >= 0);
    spqr_member member = findMember(dec, newCol->reducedMembers[id].member);
    assert(getMemberType(dec,member) == SPQR_MEMBERTYPE_RIGID);

    spqr_edge * childMarkerEdges = newCol->reducedMembers[id].childMarkerEdges;
    spqr_edge parentMarker = markerToParent(dec, member);
    spqr_node parentMarkerNodes[2] = {
            parentMarker == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeHead(dec, parentMarker),
            parentMarker == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeTail(dec, parentMarker)
    };
    spqr_node childMarkerNodes[4] = {
            childMarkerEdges[0] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeHead(dec, childMarkerEdges[0]),
            childMarkerEdges[0] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeTail(dec, childMarkerEdges[0]),
            childMarkerEdges[1] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeHead(dec, childMarkerEdges[1]),
            childMarkerEdges[1] == SPQR_INVALID_EDGE ? SPQR_INVALID_NODE : findEdgeTail(dec, childMarkerEdges[1]),
    };
    spqr_node * pathEndNodes = newCol->reducedMembers[id].rigidEndNodes;
    int numPathEndNodes = pathEndNodes[0] == SPQR_INVALID_NODE ? 0 : (pathEndNodes[2] == SPQR_INVALID_NODE ? 2 : 4);
    const int numOneEnd = newCol->reducedMembers[id].numOneEnd;
    const int numTwoEnds = newCol->reducedMembers[id].numTwoEnds;
    if(depth == 0){
        //modify the endnodes array if the parent marker is a path edge
        if(parentMarker != SPQR_INVALID_EDGE && newCol->edgeInPath[parentMarker]){
//            if(numPathEndNodes == 0) {
//                pathEndNodes[0] = parentMarkerNodes[0];
//                pathEndNodes[1] = parentMarkerNodes[1];
//                numPathEndNodes = 1;
//            }else if(numPathEndNodes == 4){
//                pathEndNodes[0] = pathEndNodes[3];
//                pathEndNodes[2] = SPQR_INVALID_NODE;
//                pathEndNodes[3] = SPQR_INVALID_NODE;
//                numPathEndNodes = 2;
//            }
        }
        assert(numPathEndNodes <= 2);
        if(numOneEnd == 0 && numTwoEnds == 0){
            //Check if an edge already exists. If there is one, we either make a parallel out of it, and the new edge
            //or if it is a parallel marker we add the edge there
            spqr_edge connectingEdge = SPQR_INVALID_EDGE;
            {
                spqr_edge iterEdge = getFirstNodeEdge(dec, pathEndNodes[0]);
                spqr_edge first = iterEdge;
                do{
                    spqr_node edgeHead = findEdgeHead(dec, iterEdge);
                    spqr_node other = edgeHead != pathEndNodes[0] ? edgeHead : findEdgeTail(dec, iterEdge);
                    if(pathEndNodes[1] == other){
                        connectingEdge = iterEdge;
                        break;
                    }
                    iterEdge = getNextNodeEdge(dec,iterEdge,pathEndNodes[0]);
                }while(iterEdge != first);
            }
            if(SPQRedgeIsInvalid(connectingEdge)){
                addTerminal(newColInfo,pathEndNodes[0]);
                addTerminal(newColInfo,pathEndNodes[1]);
                setTerminalMember(newColInfo,member);
                return MATREC_OKAY;
            }

            spqr_member parallelMember;

            spqr_member connectingMember = SPQR_INVALID_MEMBER;
            bool connectingIsParent = false;
            if(edgeIsMarker(dec,connectingEdge)){
                connectingMember = findEdgeChildMember(dec,connectingEdge);
            }else if(connectingEdge == markerToParent(dec,member)){
                connectingMember = findMemberParent(dec,member);
                assert(connectingMember != SPQR_INVALID_MEMBER);
                connectingIsParent = true;
            }

            if(SPQRmemberIsValid(connectingMember) && getMemberType(dec, connectingMember) == SPQR_MEMBERTYPE_PARALLEL){
                parallelMember = connectingMember;
            }
            else{
                //newly create the parallel member and move the edge to it
                MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &parallelMember));
                bool isTreeEdge = edgeIsTree(dec,connectingEdge);
                moveEdgeToNewMember(dec, connectingEdge, member, parallelMember);


                spqr_node oldHead = findEdgeHead(dec, connectingEdge);
                spqr_node oldTail = findEdgeTail(dec, connectingEdge);

                clearEdgeHeadAndTail(dec,connectingEdge);

                spqr_edge ignore = SPQR_INVALID_EDGE;
                spqr_edge rigidMarker = SPQR_INVALID_EDGE;
                if (connectingIsParent) {
                    MATREC_CALL(createMarkerPairWithReferences(dec, parallelMember,member, !isTreeEdge,&ignore,&rigidMarker));
                } else {
                    MATREC_CALL(createMarkerPairWithReferences(dec, member, parallelMember, isTreeEdge,&rigidMarker,&ignore));
                }
                setEdgeHeadAndTail(dec,rigidMarker,oldHead,oldTail);
            }
            setTerminalMember(newColInfo,parallelMember);
            return MATREC_OKAY;
        }
        if(numOneEnd == 1){
            spqr_node terminalNode = (pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[0] == childMarkerNodes[1])
                                     ? pathEndNodes[1] : pathEndNodes[0];
            addTerminal(newColInfo,terminalNode);

            spqr_member childMember = findEdgeChildMember(dec, childMarkerEdges[0]);
            bool headToHead = pathEndNodes[0] == childMarkerNodes[1] || pathEndNodes[1] == childMarkerNodes[1];
            assert(headToHead || pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[1] == childMarkerNodes[0]);

            MATREC_CALL(mergeMemberIntoParent(dec,childMember,headToHead));
            setTerminalMember(newColInfo, findMember(dec,member));

            return MATREC_OKAY;
        }else{
            assert(numOneEnd == 2);
            assert(newColInfo->numTerminals == 2);
            spqr_member childMember[2] = {
                    findEdgeChildMember(dec,childMarkerEdges[0]),
                    findEdgeChildMember(dec,childMarkerEdges[1])
            };
            //count to how many path end nodes the child marker is adjacent
            int numIncidentPathNodes[2] = {0,0};
            for (int c = 0; c < 2; ++c) {
                for (int i = 0; i < numPathEndNodes; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        if(pathEndNodes[i] == childMarkerNodes[2*c+j]){
                            numIncidentPathNodes[c] += 1;
                        }
                    }
                }
            }
            // If a child marker is incident to both path ends, then we ensure it is the second one.
            if(numIncidentPathNodes[0] == 2){
                spqr_member tempChildMember = childMember[0];
                childMember[0] = childMember[1];
                childMember[1] = tempChildMember;

                spqr_node tempChildMarkerNodes = childMarkerNodes[0];
                childMarkerNodes[0] = childMarkerNodes[2];
                childMarkerNodes[2] = tempChildMarkerNodes;

                tempChildMarkerNodes = childMarkerNodes[1];
                childMarkerNodes[1] = childMarkerNodes[3];
                childMarkerNodes[3] = tempChildMarkerNodes;

            }

            // Check if the two child markers are parallel. We then create a parallel with the two
            if ((childMarkerNodes[0] == childMarkerNodes[2] && childMarkerNodes[1] == childMarkerNodes[3])
                || (childMarkerNodes[0] == childMarkerNodes[3] && childMarkerNodes[1] == childMarkerNodes[2])){
                //TODO: create parallel here
                assert(false);
            }
            if(numPathEndNodes == 0){
                //We fake a path by setting the end nodes to be the common node of the child markers
                if(childMarkerNodes[0] == childMarkerNodes[2] || childMarkerNodes[0] == childMarkerNodes[3]){
                    pathEndNodes[0] = childMarkerNodes[0];
                    pathEndNodes[1] = childMarkerNodes[0];
                }else{
                    assert(childMarkerNodes[1] == childMarkerNodes[2] || childMarkerNodes[1] == childMarkerNodes[3]);
                    pathEndNodes[0] = childMarkerNodes[1];
                    pathEndNodes[1] = childMarkerNodes[1];
                }
            }
            if(pathEndNodes[0] != childMarkerNodes[0] && pathEndNodes[0] != childMarkerNodes[1]){
                spqr_node tempNode = pathEndNodes[0];
                pathEndNodes[0] = pathEndNodes[1];
                pathEndNodes[1] = tempNode;
            }
            assert(pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[0] == childMarkerNodes[1]);
            MATREC_CALL(mergeMemberIntoParent(dec,childMember[0],pathEndNodes[0] == childMarkerNodes[1]));
            MATREC_CALL(mergeMemberIntoParent(dec,childMember[1],pathEndNodes[1] == childMarkerNodes[3]));
            setTerminalMember(newColInfo, findMember(dec,member));

            return MATREC_OKAY;
        }
    }

    //non-root member

    if (numOneEnd == 0 && numTwoEnds == 0)
    {
        assert(pathEdgeIsValid(newCol->reducedMembers[id].firstPathEdge));
        assert(pathEndNodes[0] != SPQR_INVALID_NODE);
        addTerminal(newColInfo,pathEndNodes[1]);
        if (parentMarkerNodes[0] == pathEndNodes[0]){
            flipEdge(dec,parentMarker);
        }
        return MATREC_OKAY;
    }
    assert(numOneEnd == 1);

    if (numPathEndNodes >= 2) {

        if (pathEndNodes[1] != childMarkerNodes[0] && pathEndNodes[1] != childMarkerNodes[1]){
            spqr_node temp = pathEndNodes[0];
            pathEndNodes[0] = pathEndNodes[1];
            pathEndNodes[1] = temp;
        }
        assert(pathEndNodes[1] == childMarkerNodes[0] || pathEndNodes[1] == childMarkerNodes[1]);

        /* Flip parent if necessary. */
        if (pathEndNodes[0] == parentMarkerNodes[0]) {
            flipEdge(dec, parentMarker);

            spqr_node temp = parentMarkerNodes[0];
            parentMarkerNodes[0] = parentMarkerNodes[1];
            parentMarkerNodes[1] = temp;
        }

        assert(pathEndNodes[0] == parentMarkerNodes[1]);

        MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec,childMarkerEdges[0]),
                                        pathEndNodes[1] == childMarkerNodes[1]));
        return MATREC_OKAY;
    }
    /* Tested in UpdateInnerRigidNoPath. */

    /* Parent marker and child marker must be next to each other. */
    if (parentMarkerNodes[0] == childMarkerNodes[0] || parentMarkerNodes[0] == childMarkerNodes[1]){
        flipEdge(dec, parentMarker);
    }

    bool headToHead = parentMarkerNodes[0] == childMarkerNodes[1] || parentMarkerNodes[1] == childMarkerNodes[1];
    MATREC_CALL(mergeMemberIntoParent(dec, findEdgeChildMember(dec,childMarkerEdges[0]),headToHead));

    return MATREC_OKAY;
}
static MATREC_ERROR transformLoop(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition * newCol, reduced_member_id id,
                                 NewColInformation * newColInfo, int depth){
    assert(depth == 0);
    assert(dec);
    assert(newCol);
    assert(newColInfo);

    spqr_member member =newCol->reducedMembers[id].member;
    int numEdges = getNumMemberEdges(dec,member);
    assert(numEdges == 1 || numEdges == 2);
    if(getNumMemberEdges(dec,member) == 2){
        changeLoopToParallel(dec,member);
    }
    setTerminalMember(newColInfo,member);

    return MATREC_OKAY;
}
static MATREC_ERROR transformReducedMember(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition * newCol, MATRECColReducedComponent * component,
                                         reduced_member_id reducedMemberId, NewColInformation * newColInfo, int depth){
    MATRECColReducedMember * reducedMember = &newCol->reducedMembers[reducedMemberId];
    if(reducedMember->type == TYPE_CYCLE_CHILD && depth > 0){
        return MATREC_OKAY; //path has been propagated away; no need to do anything
    }

    //handle children recursively
    for (children_idx idx = reducedMember->firstChild; idx < reducedMember->firstChild + reducedMember->numChildren; ++idx){
        reduced_member_id reducedChild = newCol->childrenStorage[idx];
        assert(reducedMemberIsValid(reducedChild));
        MATREC_CALL(transformReducedMember(dec,newCol,component,reducedChild,newColInfo,depth+1));
    }
    assert(memberIsRepresentative(dec,reducedMember->member));
    SPQRMemberType type = getMemberType(dec, reducedMember->member); //TODO: find, or not?
    if(type == SPQR_MEMBERTYPE_SERIES){
        MATREC_CALL(transformSeries(dec,newCol,reducedMemberId,newColInfo,depth));
    }else if(type == SPQR_MEMBERTYPE_PARALLEL){
        MATREC_CALL(transformParallel(dec,newCol,reducedMemberId,newColInfo,depth));
    }else if(type == SPQR_MEMBERTYPE_RIGID){
        MATREC_CALL(transformRigid(dec,newCol,reducedMemberId,newColInfo,depth));
    }else if (type == SPQR_MEMBERTYPE_LOOP){
        MATREC_CALL(transformLoop(dec,newCol,reducedMemberId,newColInfo,depth));
    }
    return MATREC_OKAY;
}

static MATREC_ERROR moveReducedRoot(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol, MATRECColReducedComponent *component){
    assert(dec);
    assert(newCol);
    assert(component);

    reduced_member_id reducedMemberId = component->root;
    MATRECColReducedMember * reducedMember = &newCol->reducedMembers[reducedMemberId];
    spqr_member member = findMember(dec, reducedMember->member); //TODO: necessary find?

    SPQRMemberType type = getMemberType(dec, member);
    bool cycleWithUniqueEndChild = false;
    if(type == SPQR_MEMBERTYPE_PARALLEL){
        cycleWithUniqueEndChild = (reducedMember->numOneEnd == 1 || reducedMember->numTwoEnds == 1);
    }else if(type == SPQR_MEMBERTYPE_RIGID){
        if(reducedMember->numTwoEnds == 1 || reducedMember->numOneEnd == 1){
            spqr_node childMarkerNodes[2] = {
                    findEdgeHead(dec,reducedMember->childMarkerEdges[0]),
                    findEdgeTail(dec,reducedMember->childMarkerEdges[0]),
            };
            cycleWithUniqueEndChild = (reducedMember->rigidEndNodes[2] == SPQR_INVALID_NODE
                                       && ((reducedMember->rigidEndNodes[0] == childMarkerNodes[0]
                                            && reducedMember->rigidEndNodes[1] == childMarkerNodes[1])
                                           || (reducedMember->rigidEndNodes[0] == childMarkerNodes[1]
                                               && reducedMember->rigidEndNodes[1] == childMarkerNodes[0])));

        }else{
            cycleWithUniqueEndChild = false;
        }
    }
    if(!cycleWithUniqueEndChild){
        return MATREC_OKAY;
    }
    while(cycleWithUniqueEndChild){
        for (children_idx idx = reducedMember->firstChild; idx < reducedMember->firstChild + reducedMember->numChildren;
             ++idx) {
            reducedMemberId = newCol->childrenStorage[idx];
            ReducedMemberType childType = newCol->reducedMembers[reducedMemberId].type;
            if(childType == TYPE_SINGLE_CHILD || childType == TYPE_DOUBLE_CHILD){
                reducedMember = &newCol->reducedMembers[reducedMemberId];
                break;
            }
        }
        assert(findMemberNoCompression(dec,reducedMember->member) != member); //Assert that we found a child in the above loop

        member = reducedMember->member; //TODO: need to find, or not?
        assert(memberIsRepresentative(dec,member));

        createPathEdge(dec,newCol, markerToParent(dec,member),reducedMemberId); //Add edge which closes cycle

        assert(reducedMember->member == member);
        type = getMemberType(dec,member);
        if(type == SPQR_MEMBERTYPE_PARALLEL){
            assert(pathEdgeIsValid(reducedMember->firstPathEdge) && reducedMember->numPathEdges == 1);
            cycleWithUniqueEndChild = (reducedMember->numOneEnd == 1 || reducedMember->numTwoEnds == 1);
        }else if(type == SPQR_MEMBERTYPE_RIGID){
            if(reducedMember->numOneEnd == 1 || reducedMember->numTwoEnds == 1){
                // For non-root rigid members we have to check whether the path edges together with the parent marker edge form a cycle with a two-end child.
                spqr_edge edge = markerToParent(dec, member);
                spqr_node parentMarkerNodes[2] = {
                        findEdgeTail(dec, edge),
                        findEdgeHead(dec, edge)
                };
                spqr_node childMarkerNodes[2] = {
                        findEdgeTail(dec, reducedMember->childMarkerEdges[0]),
                        findEdgeHead(dec, reducedMember->childMarkerEdges[0])
                };

                int numEndNodes = reducedMember->rigidEndNodes[0] == SPQR_INVALID_NODE ?
                                  0 : (reducedMember->rigidEndNodes[2] == SPQR_INVALID_NODE ? 2 : 4);
                if (numEndNodes == 0)
                {
                    /* Without path edges the child marker would have to be parallel to the parent marker, which is detected
                     * during typing. */
                    cycleWithUniqueEndChild = false;
                    break;
                }

                /* Determine the end nodes of the path including the parent marker edge. */
                spqr_node endNodes[2];
                if (numEndNodes == 4)
                {
                    endNodes[0] = reducedMember->rigidEndNodes[1];
                    endNodes[1] = reducedMember->rigidEndNodes[3];
                }
                else if (reducedMember->rigidEndNodes[0] == parentMarkerNodes[0])
                {
                    endNodes[0] = parentMarkerNodes[1];
                    endNodes[1] = reducedMember->rigidEndNodes[1];
                }
                else
                {
                    assert(reducedMember->rigidEndNodes[0] == parentMarkerNodes[1]);
                    endNodes[0] = parentMarkerNodes[0];
                    endNodes[1] = reducedMember->rigidEndNodes[1];
                }

                cycleWithUniqueEndChild = (endNodes[0] == childMarkerNodes[0] && endNodes[1] == childMarkerNodes[1])
                                          || (endNodes[0] == childMarkerNodes[1] && endNodes[1] == childMarkerNodes[0]);
            }else{
                cycleWithUniqueEndChild = false;
            }
        }else{
            assert(type == SPQR_MEMBERTYPE_SERIES);
            if(reducedMember->numOneEnd == 1 || reducedMember->numTwoEnds == 1){
                cycleWithUniqueEndChild = reducedMember->numPathEdges == (getNumMemberEdges(dec,member) -1);
            }else{
                cycleWithUniqueEndChild = false;
            }
        }
    }
    component->root = reducedMemberId;
    return MATREC_OKAY;
}
static MATREC_ERROR transformComponent(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol, MATRECColReducedComponent * component, NewColInformation * newColInfo){
    assert(dec);
    assert(newCol);
    assert(component);
    assert(newColInfo);

    //First, ensure that the reduced root member is moved up if it is a cycle
    MATREC_CALL(moveReducedRoot(dec,newCol,component));
    //Then, recursively transform the components
    MATREC_CALL(transformReducedMember(dec,newCol,component,component->root,newColInfo,0));
    return MATREC_OKAY;
}

MATREC_ERROR MATRECGraphicColumnAdditionAdd(MATRECGraphicDecomposition *dec, MATRECGraphicColumnAddition *newCol){
    assert(dec);
    assert(newCol);

    if(newCol->numReducedComponents == 0){
        spqr_member member;
        MATREC_CALL(createStandaloneSeries(dec,newCol->newRowEdges,newCol->numNewRowEdges,newCol->newColIndex,&member));
    }else if(newCol->numReducedComponents == 1){
        NewColInformation information = emptyNewColInformation();
        MATREC_CALL(transformComponent(dec,newCol,&newCol->reducedComponents[0],&information));
        assert(memberIsRepresentative(dec,information.member));
        if(newCol->numNewRowEdges == 0){
            spqr_edge colEdge = SPQR_INVALID_EDGE;
            MATREC_CALL(createColumnEdge(dec,information.member,&colEdge,newCol->newColIndex));
            if(SPQRnodeIsValid(information.terminalNode[0])){
                setEdgeHeadAndTail(dec,colEdge,
                                   findNode(dec,information.terminalNode[0]),findNode(dec,information.terminalNode[1]));
            }
        }else{
            spqr_member newSeries;
            MATREC_CALL(createConnectedSeries(dec,newCol->newRowEdges,newCol->numNewRowEdges,newCol->newColIndex,&newSeries));
            spqr_edge markerEdge = SPQR_INVALID_EDGE;
            spqr_edge ignore = SPQR_INVALID_EDGE;
            MATREC_CALL(createMarkerPairWithReferences(dec,information.member,newSeries,false,&markerEdge,&ignore));
            if(SPQRnodeIsValid(information.terminalNode[0])){
                setEdgeHeadAndTail(dec,markerEdge,findNode(dec,information.terminalNode[0]),
                                   findNode(dec,information.terminalNode[1]));
            }
        }
        if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
            assert(getNumMemberEdges(dec,information.member) == 2 || getNumMemberEdges(dec,information.member) == 3);
            if(getNumMemberEdges(dec,information.member) == 3){
                changeLoopToParallel(dec,information.member);
            }
        }
    }else{
#ifndef NDEBUG
        int numDecComponentsBefore = numConnectedComponents(dec);
#endif
        spqr_member newSeries;
        MATREC_CALL(createConnectedSeries(dec,newCol->newRowEdges,newCol->numNewRowEdges,newCol->newColIndex,&newSeries));
        for (int i = 0; i < newCol->numReducedComponents; ++i) {
            NewColInformation information = emptyNewColInformation();
            MATREC_CALL(transformComponent(dec,newCol,&newCol->reducedComponents[i],&information));
            if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
                assert(getNumMemberEdges(dec,information.member) == 1);
                moveEdgeToNewMember(dec, getFirstMemberEdge(dec,information.member),information.member,newSeries);
                dec->members[information.member].type = SPQR_MEMBERTYPE_UNASSIGNED;
            }else {
                reorderComponent(dec,information.member); //reorder the subtree so that the new series member is a parent
                spqr_edge markerEdge = SPQR_INVALID_EDGE;
                spqr_edge ignore = SPQR_INVALID_EDGE;
                MATREC_CALL(
                        createMarkerPairWithReferences(dec, newSeries, information.member, true, &ignore, &markerEdge));
                if (SPQRnodeIsValid(information.terminalNode[0])) {
                    setEdgeHeadAndTail(dec, markerEdge, findNode(dec, information.terminalNode[0]),
                                       findNode(dec, information.terminalNode[1]));
                }
            }
        }
        decreaseNumConnectedComponents(dec,newCol->numReducedComponents-1);
        assert(numConnectedComponents(dec) == (numDecComponentsBefore - newCol->numReducedComponents + 1));
    }
    return MATREC_OKAY;
}

bool MATRECGraphicColumnAdditionRemainsGraphic(MATRECGraphicColumnAddition *newCol){
    return newCol->remainsGraphic;
}


static int min(int a, int b){
    return a < b ? a : b;
}

typedef int cut_edge_id;
#define INVALID_CUT_EDGE (-1)

static bool cutEdgeIsInvalid(const cut_edge_id edge){
    return edge < 0;
}
static bool cutEdgeIsValid(const cut_edge_id edge){
    return !cutEdgeIsInvalid(edge);
}

typedef struct { //TODO:test if overhead of pointers is worth it?
    spqr_edge edge;
    cut_edge_id nextMember;
    cut_edge_id nextOverall;
} CutEdgeListNode;

typedef enum{
    TYPE_UNDETERMINED = 0,
    TYPE_PROPAGATED = 1,
    TYPE_MERGED = 2,
    TYPE_NOT_GRAPHIC = 3
} RowReducedMemberType;

/**
 * A 'set' of two nodes. We need this a lot in the used algorithms
 */
typedef struct {
    spqr_node first;
    spqr_node second;
}NodePair;

static void NodePairEmptyInitialize(NodePair * pair){
    pair->first = SPQR_INVALID_NODE;
    pair->second = SPQR_INVALID_NODE;
}
static void NodePairInitialize(NodePair * pair, const spqr_node firstNode, const spqr_node secondNode){
    pair->first = firstNode;
    pair->second = secondNode;
}
static bool NodePairIsEmpty(const NodePair * pair){
    return SPQRnodeIsInvalid(pair->first);
}
static bool NodePairHasTwo(const NodePair * pair){
    return SPQRnodeIsValid(pair->second);
}
static void NodePairIntersection(NodePair * toChange, const spqr_node first, const spqr_node second){
    if (toChange->first!= first &&  toChange->first != second) {
        toChange->first = SPQR_INVALID_NODE;
    }
    if (toChange->second != first && toChange->second != second) {
        toChange->second= SPQR_INVALID_NODE;
    }
    if (SPQRnodeIsInvalid(toChange->first) && SPQRnodeIsValid(toChange->second)) {
        swap_ints(&toChange->first,&toChange->second);
    }
}
static void NodePairInsert(NodePair * pair, const spqr_node node){
    if(pair->first != SPQR_INVALID_NODE){
        assert(pair->second == SPQR_INVALID_NODE);
        pair->second = node;
    }else{
        pair->first = node;
    }
}

typedef struct{
    int low;
    int discoveryTime;
} ArticulationNodeInformation;

//We allocate the callstacks of recursive algorithms (usually DFS, bounded by some linear number of calls)
//If one does not do this, we overflow the stack for large matrices/graphs through the number of recursive function calls
//Then, we can write the recursive algorithms as while loops and allocate the function call data on the heap, preventing
//Stack overflows
typedef struct {
    spqr_node node;
    spqr_edge nodeEdge;
} DFSCallData;

typedef struct{
    children_idx currentChild;
    reduced_member_id id;
} MergeTreeCallData;

typedef struct{
    spqr_node node;
    spqr_edge edge;
} ColorDFSCallData;

typedef struct{
    spqr_edge edge;
    spqr_node node;
    spqr_node parent;
    bool isAP;
} ArticulationPointCallStack;

typedef enum{
    UNCOLORED = 0,
    COLOR_FIRST = 1,
    COLOR_SECOND = 2
} COLOR_STATUS;

typedef struct {
    spqr_member member;
    spqr_member rootMember;
    int depth;
    RowReducedMemberType type;
    reduced_member_id parent;

    children_idx firstChild;
    children_idx numChildren;
    children_idx numPropagatedChildren;

    cut_edge_id firstCutEdge;
    int numCutEdges;

    NodePair splitting_nodes;
    bool allHaveCommonNode;
    spqr_edge articulationEdge;

    spqr_node coloredNode; //points to a colored node so that we can efficiently zero out the colors again.
} MATRECRowReducedMember;

typedef struct {
    int rootDepth;
    reduced_member_id root;
} MATRECRowReducedComponent;

struct MATRECGraphicRowAdditionImpl{
    bool remainsGraphic;

    MATRECRowReducedMember* reducedMembers;
    int memReducedMembers;
    int numReducedMembers;

    MATRECRowReducedComponent* reducedComponents;
    int memReducedComponents;
    int numReducedComponents;

    MemberInfo * memberInformation;
    int memMemberInformation;
    int numMemberInformation;

    reduced_member_id * childrenStorage;
    int memChildrenStorage;
    int numChildrenStorage;

    CutEdgeListNode * cutEdges;
    int memCutEdges;
    int numCutEdges;
    cut_edge_id firstOverallCutEdge;

    MATREC_row newRowIndex;

    MATREC_col * newColumnEdges;
    int memColumnEdges;
    int numColumnEdges;

    reduced_member_id * leafMembers;
    int numLeafMembers;
    int memLeafMembers;

    spqr_edge * decompositionColumnEdges;
    int memDecompositionColumnEdges;
    int numDecompositionColumnEdges;

    bool * isEdgeCut;
    int numIsEdgeCut;
    int memIsEdgeCut;

    COLOR_STATUS * nodeColors;
    int memNodeColors;

    spqr_node * articulationNodes;
    int numArticulationNodes;
    int memArticulationNodes;

    ArticulationNodeInformation * articulationNodeSearchInfo;
    int memNodeSearchInfo;

    int * crossingPathCount;
    int memCrossingPathCount;

    DFSCallData * intersectionDFSData;
    int memIntersectionDFSData;

    ColorDFSCallData * colorDFSData;
    int memColorDFSData;

    ArticulationPointCallStack * artDFSData;
    int memArtDFSData;

    CreateReducedMembersCallstack * createReducedMembersCallstack;
    int memCreateReducedMembersCallstack;

    int * intersectionPathDepth;
    int memIntersectionPathDepth;

    spqr_node * intersectionPathParent;
    int memIntersectionPathParent;

    MergeTreeCallData * mergeTreeCallData;
    int memMergeTreeCallData;

};

typedef struct {
    spqr_member member;
    spqr_node firstNode;
    spqr_node secondNode;
} NewRowInformation;

static NewRowInformation emptyNewRowInformation(void){
    NewRowInformation information;
    information.member = SPQR_INVALID_MEMBER;
    information.firstNode = SPQR_INVALID_NODE;
    information.secondNode = SPQR_INVALID_NODE;
    return information;
}

/**
 * Saves the information of the current row and partitions it based on whether or not the given columns are
 * already part of the decomposition.
 */
static MATREC_ERROR newRowUpdateRowInformation(const MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow, const MATREC_row row, const MATREC_col * columns, const size_t numColumns)
{
    newRow->newRowIndex = row;

    newRow->numDecompositionColumnEdges = 0;
    newRow->numColumnEdges = 0;

    for (size_t i = 0; i < numColumns; ++i) {
        spqr_edge columnEdge = getDecompositionColumnEdge(dec, columns[i]);
        if(SPQRedgeIsValid(columnEdge)){ //If the edge is the current decomposition: save it in the array
            if(newRow->numDecompositionColumnEdges == newRow->memDecompositionColumnEdges){
                int newNumEdges = newRow->memDecompositionColumnEdges == 0 ? 8 : 2*newRow->memDecompositionColumnEdges; //TODO: make reallocation numbers more consistent with rest?
                newRow->memDecompositionColumnEdges = newNumEdges;
                MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->decompositionColumnEdges,
                                                (size_t) newRow->memDecompositionColumnEdges));
            }
            newRow->decompositionColumnEdges[newRow->numDecompositionColumnEdges] = columnEdge;
            ++newRow->numDecompositionColumnEdges;
        }else{
            //Not in the decomposition: add it to the set of edges which are newly added with this row.
            if(newRow->numColumnEdges == newRow->memColumnEdges){
                int newNumEdges = newRow->memColumnEdges == 0 ? 8 : 2*newRow->memColumnEdges; //TODO: make reallocation numbers more consistent with rest?
                newRow->memColumnEdges = newNumEdges;
                MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->newColumnEdges,
                                                (size_t)newRow->memColumnEdges));
            }
            newRow->newColumnEdges[newRow->numColumnEdges] = columns[i];
            newRow->numColumnEdges++;
        }
    }

    return MATREC_OKAY;
}

/**
 * Recursively creates reduced members from this member to the root of the decomposition tree.
 * @param dec
 * @param newRow
 * @param member
 * @return
 */
static reduced_member_id createRowReducedMembersToRoot(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition * newRow, const spqr_member firstMember){
    assert(SPQRmemberIsValid(firstMember));

    CreateReducedMembersCallstack * callstack = newRow->createReducedMembersCallstack;
    callstack[0].member = firstMember;
    int callDepth = 0;

    while(callDepth >= 0){
        spqr_member member = callstack[callDepth].member;
        reduced_member_id reducedMember = newRow->memberInformation[member].reducedMember;

        bool reducedValid = reducedMemberIsValid(reducedMember);
        if(!reducedValid) {
            //reduced member was not yet created; we create it
            reducedMember = newRow->numReducedMembers;

            MATRECRowReducedMember *reducedMemberData = &newRow->reducedMembers[reducedMember];
            ++newRow->numReducedMembers;

            reducedMemberData->member = member;
            reducedMemberData->numChildren = 0;
            reducedMemberData->numCutEdges = 0;
            reducedMemberData->firstCutEdge = INVALID_CUT_EDGE;
            reducedMemberData->type = TYPE_UNDETERMINED;
            reducedMemberData->numPropagatedChildren = 0;
            reducedMemberData->allHaveCommonNode = false;
            reducedMemberData->articulationEdge = SPQR_INVALID_EDGE;
            reducedMemberData->coloredNode = SPQR_INVALID_NODE;
            NodePairEmptyInitialize(&reducedMemberData->splitting_nodes);

            newRow->memberInformation[member].reducedMember = reducedMember;
            assert(memberIsRepresentative(dec, member));
            spqr_member parentMember = findMemberParent(dec, member);

            if (SPQRmemberIsValid(parentMember)) {
                //recursive call to parent member
                ++callDepth;
                assert(callDepth < newRow->memCreateReducedMembersCallstack);
                callstack[callDepth].member = parentMember;
                continue;

            } else {
                //we found a new reduced decomposition component

                reducedMemberData->parent = INVALID_REDUCED_MEMBER;
                reducedMemberData->depth = 0;
                reducedMemberData->rootMember = member;

                assert(newRow->numReducedComponents < newRow->memReducedComponents);
                newRow->reducedComponents[newRow->numReducedComponents].root = reducedMember;
                ++newRow->numReducedComponents;
            }
        }
        if(reducedValid){
            assert(reducedMember < newRow->numReducedMembers);
            //Reduced member was already created in earlier call
            //update the depth of the root if appropriate
            reduced_member_id * depthMinimizer = &newRow->memberInformation[newRow->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
            if(reducedMemberIsInvalid(*depthMinimizer) ||
               newRow->reducedMembers[reducedMember].depth < newRow->reducedMembers[*depthMinimizer].depth){
                *depthMinimizer = reducedMember;
            }
        }
        while(true){
            --callDepth;
            if(callDepth < 0 ) break;
            spqr_member parentMember = callstack[callDepth + 1].member;
            reduced_member_id parentReducedMember = newRow->memberInformation[parentMember].reducedMember;
            spqr_member currentMember = callstack[callDepth].member;
            reduced_member_id currentReducedMember = newRow->memberInformation[currentMember].reducedMember;

            MATRECRowReducedMember *parentReducedMemberData = &newRow->reducedMembers[parentReducedMember];
            MATRECRowReducedMember *reducedMemberData = &newRow->reducedMembers[currentReducedMember];

            reducedMemberData->parent = parentReducedMember;
            reducedMemberData->depth = parentReducedMemberData->depth + 1;
            reducedMemberData->rootMember = parentReducedMemberData->rootMember;

            newRow->reducedMembers[parentReducedMember].numChildren++;
        }

    }

    reduced_member_id returnedMember = newRow->memberInformation[callstack[0].member].reducedMember;
    return returnedMember;
}


/**
 * Construct a smaller sub tree of the decomposition on which the cut edges lie.
 * @return
 */
static MATREC_ERROR constructRowReducedDecomposition(MATRECGraphicDecomposition* dec, MATRECGraphicRowAddition* newRow){
    //TODO: chop up into more functions
    //TODO: stricter assertions/array bounds checking in this function
#ifndef NDEBUG
    for (int i = 0; i < newRow->memMemberInformation; ++i) {
        assert(reducedMemberIsInvalid(newRow->memberInformation[i].reducedMember));
    }
#endif

    newRow->numReducedComponents = 0;
    newRow->numReducedMembers = 0;
    if(newRow->numDecompositionColumnEdges == 0){ //Early return in case the reduced decomposition will be empty
        return MATREC_OKAY;
    }
    assert(newRow->numReducedMembers == 0);
    assert(newRow->numReducedComponents == 0);

    int newSize = largestMemberID(dec); //Is this sufficient?
    if(newSize > newRow->memReducedMembers){
        newRow->memReducedMembers = max(2*newRow->memReducedMembers,newSize);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->reducedMembers,(size_t) newRow->memReducedMembers));
    }
    if(newSize > newRow->memMemberInformation){
        int updatedSize = max(2*newRow->memMemberInformation,newSize);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->memberInformation,(size_t) updatedSize));
        for (int i = newRow->memMemberInformation; i < updatedSize; ++i) {
            newRow->memberInformation[i].reducedMember = INVALID_REDUCED_MEMBER;
            newRow->memberInformation[i].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
        }
        newRow->memMemberInformation = updatedSize;

    }

    int numComponents = numConnectedComponents(dec);
    if(numComponents > newRow->memReducedComponents){
        newRow->memReducedComponents = max(2*newRow->memReducedComponents,numComponents);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->reducedComponents,(size_t) newRow->memReducedComponents));
    }

    int numMembers = getNumMembers(dec);
    if(newRow->memCreateReducedMembersCallstack < numMembers){
        newRow->memCreateReducedMembersCallstack = max(2*newRow->memCreateReducedMembersCallstack,numMembers);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->createReducedMembersCallstack,(size_t) newRow->memCreateReducedMembersCallstack));
    }

    //Create the reduced members (recursively)
    for (int i = 0; i < newRow->numDecompositionColumnEdges; ++i) {
        assert(i < newRow->memDecompositionColumnEdges);
        spqr_edge edge = newRow->decompositionColumnEdges[i];
        spqr_member edgeMember = findEdgeMember(dec, edge);
        reduced_member_id reducedMember = createRowReducedMembersToRoot(dec,newRow,edgeMember);
        reduced_member_id* depthMinimizer = &newRow->memberInformation[newRow->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
        if(reducedMemberIsInvalid(*depthMinimizer)){
            *depthMinimizer = reducedMember;
        }
    }

    //Set the reduced roots according to the root depth minimizers
    for (int i = 0; i < newRow->numReducedComponents; ++i) {
        MATRECRowReducedComponent * component = &newRow->reducedComponents[i];
        spqr_member rootMember = newRow->reducedMembers[component->root].member;
        reduced_member_id reducedMinimizer = newRow->memberInformation[rootMember].rootDepthMinimizer;
        component->rootDepth =  newRow->reducedMembers[reducedMinimizer].depth;
        component->root = reducedMinimizer;

        //This simplifies code further down which does not need to be component-aware; just pretend that the reduced member is the new root.
        newRow->reducedMembers[component->root].parent = INVALID_REDUCED_MEMBER;
        assert(memberIsRepresentative(dec,rootMember));
    }

    //update the children array
    int numTotalChildren = 0;
    for (int i = 0; i < newRow->numReducedMembers; ++i) {
        MATRECRowReducedMember * reducedMember = &newRow->reducedMembers[i];
        reduced_member_id minimizer = newRow->memberInformation[reducedMember->rootMember].rootDepthMinimizer;
        if(reducedMember->depth >= newRow->reducedMembers[minimizer].depth){
            reducedMember->firstChild = numTotalChildren;
            numTotalChildren += reducedMember->numChildren;
            reducedMember->numChildren = 0;
        }
    }

    if(newRow->memChildrenStorage < numTotalChildren){
        int newMemSize = max(newRow->memChildrenStorage*2, numTotalChildren);
        newRow->memChildrenStorage = newMemSize;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->childrenStorage,(size_t) newRow->memChildrenStorage));
    }
    newRow->numChildrenStorage = numTotalChildren;

    //Fill up the children array`
    for(reduced_member_id  reducedMember = 0; reducedMember < newRow->numReducedMembers; ++reducedMember){
        MATRECRowReducedMember * reducedMemberData = &newRow->reducedMembers[reducedMember];
        if(reducedMemberData->depth <= newRow->reducedMembers[newRow->memberInformation[reducedMemberData->rootMember].rootDepthMinimizer].depth){
            continue;
        }
        spqr_member parentMember = findMemberParent(dec, reducedMemberData->member);
        reduced_member_id parentReducedMember = SPQRmemberIsValid(parentMember) ? newRow->memberInformation[parentMember].reducedMember : INVALID_REDUCED_MEMBER;
        if(reducedMemberIsValid(parentReducedMember)){ //TODO: probably one of these two checks/branches is unnecessary, as there is a single failure case? (Not sure)
            MATRECRowReducedMember * parentReducedMemberData = &newRow->reducedMembers[parentReducedMember];
            newRow->childrenStorage[parentReducedMemberData->firstChild + parentReducedMemberData->numChildren] = reducedMember;
            ++parentReducedMemberData->numChildren;
        }
    }

    //Clean up the root depth minimizers.
    for (int i = 0; i < newRow->numReducedMembers; ++i) {
        MATRECRowReducedMember * reducedMember = &newRow->reducedMembers[i];
        assert(reducedMember);
        spqr_member rootMember = reducedMember->rootMember;
        assert(rootMember >= 0);
        assert(rootMember < dec->memMembers);
        newRow->memberInformation[rootMember].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
    }

    return MATREC_OKAY;
}


/**
 * Marks an edge as 'cut'. This implies that its cycle in the decomposition must be elongated
 * @param newRow
 * @param edge
 * @param reducedMember
 */
static void createCutEdge(MATRECGraphicRowAddition* newRow, const spqr_edge edge, const reduced_member_id reducedMember){
    cut_edge_id cut_edge =  newRow->numCutEdges;

    CutEdgeListNode * listNode = &newRow->cutEdges[cut_edge];
    listNode->edge = edge;

    listNode->nextMember = newRow->reducedMembers[reducedMember].firstCutEdge;
    newRow->reducedMembers[reducedMember].firstCutEdge = cut_edge;

    listNode->nextOverall = newRow->firstOverallCutEdge;
    newRow->firstOverallCutEdge = cut_edge;

    newRow->numCutEdges++;
    newRow->reducedMembers[reducedMember].numCutEdges++;
    assert(newRow->numCutEdges <= newRow->memCutEdges);

    assert(edge < newRow->memIsEdgeCut);
    newRow->isEdgeCut[edge] = true;
}

/**
 * Creates all cut edges within the decomposition for the new row.
 * Note this preallocates memory for cut edges which may be created by propagation.
 */
static MATREC_ERROR createReducedDecompositionCutEdges(MATRECGraphicDecomposition* dec, MATRECGraphicRowAddition* newRow){
    //Allocate memory for cut edges
    spqr_edge maxEdgeID = largestEdgeID(dec);
    if(maxEdgeID > newRow->memIsEdgeCut){
        int newSize = max(maxEdgeID,2*newRow->memIsEdgeCut);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->isEdgeCut,(size_t) newSize));
        for (int i = newRow->memIsEdgeCut; i < newSize ; ++i) {
            newRow->isEdgeCut[i] = false;
        }
        newRow->memIsEdgeCut = newSize;
    }
#ifndef NDEBUG
    for (int i = 0; i < newRow->memIsEdgeCut; ++i) {
        assert(!newRow->isEdgeCut[i]);
    }
#endif

    int numNeededEdges = newRow->numDecompositionColumnEdges*4; //3 Is not enough; see tests. Probably 3 + 12 or so is, but cannot be bothered to work that out for now
    if(numNeededEdges > newRow->memCutEdges){
        int newSize = max(newRow->memCutEdges*2, numNeededEdges);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->cutEdges,(size_t) newSize));
        newRow->memCutEdges = newSize;
    }
    newRow->numCutEdges = 0;
    newRow->firstOverallCutEdge = INVALID_CUT_EDGE;
    for (int i = 0; i < newRow->numDecompositionColumnEdges; ++i) {
        spqr_edge edge = newRow->decompositionColumnEdges[i];
        spqr_member member = findEdgeMember(dec, edge);
        reduced_member_id reduced_member = newRow->memberInformation[member].reducedMember;
        assert(reducedMemberIsValid(reduced_member));
        createCutEdge(newRow,edge,reduced_member);
    }

    return MATREC_OKAY;
}

/**
 * Determines the members of the reduced decomposition which are leafs.
 * This is used in propagation to ensure propagation is only checked for components which have at most one neighbour
 * which is not propagated.
 */
static MATREC_ERROR determineLeafReducedMembers(const MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow) {
    if (newRow->numDecompositionColumnEdges > newRow->memLeafMembers) {
        newRow->memLeafMembers = max(newRow->numDecompositionColumnEdges, 2 * newRow->memLeafMembers);
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &newRow->leafMembers, (size_t) newRow->memLeafMembers));
    }
    newRow->numLeafMembers = 0;

    for (reduced_member_id reducedMember = 0; reducedMember < newRow->numReducedMembers; ++reducedMember) {
        if (newRow->reducedMembers[reducedMember].numChildren == 0) {
            newRow->leafMembers[newRow->numLeafMembers] = reducedMember;
            ++newRow->numLeafMembers;
        }
    }
    return MATREC_OKAY;
}

/**
 * Preallocates memory arrays necessary for searching rigid components.
 */
static MATREC_ERROR allocateRigidSearchMemory(const MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow){
    int totalNumNodes = getNumNodes(dec);
    if(totalNumNodes > newRow->memNodeColors){
        int newSize = max(2*newRow->memNodeColors,totalNumNodes);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->nodeColors,(size_t) newSize));
        for (int i = newRow->memNodeColors; i < newSize; ++i) {
            newRow->nodeColors[i] = UNCOLORED;
        }
        newRow->memNodeColors = newSize;
    }

    if(totalNumNodes > newRow->memArticulationNodes){
        int newSize = max(2*newRow->memArticulationNodes,totalNumNodes);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->articulationNodes,(size_t) newSize));
        newRow->memArticulationNodes = newSize;
    }
    if(totalNumNodes > newRow->memNodeSearchInfo){
        int newSize = max(2*newRow->memNodeSearchInfo,totalNumNodes);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->articulationNodeSearchInfo,(size_t) newSize));
        newRow->memNodeSearchInfo = newSize;
    }
    if(totalNumNodes > newRow->memCrossingPathCount){
        int newSize = max(2*newRow->memCrossingPathCount,totalNumNodes);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->crossingPathCount,(size_t) newSize));
        newRow->memCrossingPathCount = newSize;
    }

    //TODO: see if tradeoff for performance bound by checking max # of nodes of rigid is worth it to reduce size
    //of the following allocations
    int largestID  = largestNodeID(dec); //TODO: only update the stack sizes of the following when needed? The preallocation might be causing performance problems
    if(largestID > newRow->memIntersectionDFSData){
        int newSize = max(2*newRow->memIntersectionDFSData,largestID);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->intersectionDFSData,(size_t) newSize));
        newRow->memIntersectionDFSData = newSize;
    }
    if(largestID > newRow->memColorDFSData){
        int newSize = max(2*newRow->memColorDFSData, largestID);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->colorDFSData,(size_t) newSize));
        newRow->memColorDFSData = newSize;
    }
    if(largestID > newRow->memArtDFSData){
        int newSize = max(2*newRow->memArtDFSData,largestID);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->artDFSData,(size_t) newSize));
        newRow->memArtDFSData = newSize;
    }

    for (int i = 0; i < newRow->memIntersectionPathDepth; ++i) {
        newRow->intersectionPathDepth[i] = -1;
    }

    if(largestID > newRow->memIntersectionPathDepth){
        int newSize = max(2*newRow->memIntersectionPathDepth,largestID);
        MATRECreallocBlockArray(dec->env, &newRow->intersectionPathDepth, (size_t) newSize);
        for (int i = newRow->memIntersectionPathDepth; i < newSize; ++i) {
            newRow->intersectionPathDepth[i] = -1;
        }
        newRow->memIntersectionPathDepth = newSize;
    }
    for (int i = 0; i < newRow->memIntersectionPathParent; ++i) {
        newRow->intersectionPathParent[i] = SPQR_INVALID_NODE;
    }
    if(largestID > newRow->memIntersectionPathParent){
        int newSize = max(2*newRow->memIntersectionPathParent,largestID);
        MATRECreallocBlockArray(dec->env, &newRow->intersectionPathParent, (size_t) newSize);
        for (int i = newRow->memIntersectionPathParent; i <newSize; ++i) {
            newRow->intersectionPathParent[i] = SPQR_INVALID_NODE;
        }
        newRow->memIntersectionPathParent = newSize;
    }

    return MATREC_OKAY;
}
static void zeroOutColors(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow, const spqr_node firstRemoveNode){
    assert(firstRemoveNode < newRow->memNodeColors);

    newRow->nodeColors[firstRemoveNode] = UNCOLORED;
    ColorDFSCallData * data = newRow->colorDFSData;

    data[0].node = firstRemoveNode;
    data[0].edge = getFirstNodeEdge(dec,firstRemoveNode);

    int depth = 0;

    while(depth >= 0){
        assert(depth < newRow->memColorDFSData);
        ColorDFSCallData * callData = &data[depth];
        spqr_node head = findEdgeHead(dec, callData->edge);
        spqr_node tail = findEdgeTail(dec, callData->edge);
        spqr_node otherNode = callData->node == head ? tail : head;
        assert(otherNode < newRow->memNodeColors);
        if(newRow->nodeColors[otherNode] != UNCOLORED){
            callData->edge = getNextNodeEdge(dec,callData->edge,callData->node);

            newRow->nodeColors[otherNode] = UNCOLORED;
            ++depth;
            data[depth].node = otherNode;
            data[depth].edge = getFirstNodeEdge(dec,otherNode);
            continue;
        }

        callData->edge = getNextNodeEdge(dec,callData->edge,callData->node);
        while(depth >= 0 && data[depth].edge == getFirstNodeEdge(dec,data[depth].node)){
            --depth;
        }
    }


}

static void cleanUpPreviousIteration(MATRECGraphicDecomposition * dec, MATRECGraphicRowAddition * newRow){
    //zero out coloring information from previous check
    for (int i = 0; i < newRow->numReducedMembers; ++i) {
        if(SPQRnodeIsValid(newRow->reducedMembers[i].coloredNode)){
            zeroOutColors(dec,newRow,newRow->reducedMembers[i].coloredNode);
            newRow->reducedMembers[i].coloredNode = SPQR_INVALID_NODE;
        }
    }
#ifndef NDEBUG
    for (int i = 0; i < newRow->memNodeColors; ++i) {
        assert(newRow->nodeColors[i] == UNCOLORED);
    }
#endif

    //For cut edges: clear them from the array from previous iteration
    cut_edge_id cutEdgeIdx = newRow->firstOverallCutEdge;
    while(cutEdgeIsValid(cutEdgeIdx)){
        spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
        cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextOverall;
        newRow->isEdgeCut[cutEdge] = false;
    }
}

static NodePair rigidDetermineAllAdjacentSplittableNodes(MATRECGraphicDecomposition * dec, MATRECGraphicRowAddition * newRow,
                                                         const reduced_member_id toCheck, spqr_edge * secondArticulationEdge){
    NodePair pair;
    NodePairEmptyInitialize(&pair);
    NodePair * nodePair = &pair;
    assert(NodePairIsEmpty(nodePair));
    assert(newRow->reducedMembers[toCheck].numCutEdges > 0);//calling this function otherwise is nonsensical
    assert(*secondArticulationEdge == SPQR_INVALID_EDGE);

    cut_edge_id cutEdgeIdx = newRow->reducedMembers[toCheck].firstCutEdge;
    spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
    spqr_node head = findEdgeHead(dec, cutEdge);
    spqr_node tail = findEdgeTail(dec, cutEdge);
    NodePairInitialize(nodePair, head, tail);

    while (cutEdgeIsValid(newRow->cutEdges[cutEdgeIdx].nextMember)) {
        cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextMember;
        cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
        head = findEdgeHead(dec, cutEdge);
        tail = findEdgeTail(dec, cutEdge);
        NodePairIntersection(nodePair, head, tail);

        if (NodePairIsEmpty(nodePair)) {
            break;
        }
    }

    if(!NodePairIsEmpty(nodePair)){
        //Check if the cut edges are n-1 of the n edges incident at this node; if so, that point is an articulation point
        //as it disconnects the first splitnode from the rest of the tree
        spqr_node splitNode = nodePair->first;
        if(!NodePairHasTwo(nodePair) &&newRow->reducedMembers[toCheck].numCutEdges == nodeDegree(dec,splitNode) -1 ){
            spqr_edge firstNodeEdge = getFirstNodeEdge(dec, splitNode);
            spqr_edge neighbourEdge = firstNodeEdge;
            do{
                if(edgeIsTree(dec,neighbourEdge)){
                    break;
                }
                neighbourEdge = getNextNodeEdge(dec,neighbourEdge,splitNode);
            }while(neighbourEdge != firstNodeEdge);
            spqr_node otherHead = findEdgeHead(dec, neighbourEdge);
            spqr_node otherTail = findEdgeTail(dec, neighbourEdge);
            spqr_node otherNode = otherHead == splitNode ? otherTail : otherHead;
            NodePairInsert(nodePair,otherNode);
            *secondArticulationEdge = neighbourEdge;
        }
    }
    return pair;
}

//TODO: remove MATREC_ERROR from below functions (until propagation function, basically) and refactor memory allocation
static MATREC_ERROR zeroOutColorsExceptNeighbourhood(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                                   const spqr_node articulationNode, const spqr_node startRemoveNode){
    COLOR_STATUS * neighbourColors;
    int degree = nodeDegree(dec,articulationNode);
    MATREC_CALL(MATRECallocBlockArray(dec->env,&neighbourColors,(size_t) degree));

    {
        int i = 0;
        spqr_edge artFirstEdge = getFirstNodeEdge(dec, articulationNode);
        spqr_edge artItEdge = artFirstEdge;
        do{
            spqr_node head = findEdgeHead(dec, artItEdge);
            spqr_node tail = findEdgeTail(dec, artItEdge);
            spqr_node otherNode = articulationNode == head ? tail : head;
            neighbourColors[i] = newRow->nodeColors[otherNode];
            i++;
            assert(i <= degree);
            artItEdge = getNextNodeEdge(dec,artItEdge,articulationNode);
        }while(artItEdge != artFirstEdge);
    }
    zeroOutColors(dec,newRow,startRemoveNode);

    {
        int i = 0;
        spqr_edge artFirstEdge = getFirstNodeEdge(dec, articulationNode);
        spqr_edge artItEdge = artFirstEdge;
        do{
            spqr_node head = findEdgeHead(dec, artItEdge);
            spqr_node tail = findEdgeTail(dec, artItEdge);
            spqr_node otherNode = articulationNode == head ? tail : head;
            newRow->nodeColors[otherNode] = neighbourColors[i];
            i++;
            assert(i <= degree);
            artItEdge = getNextNodeEdge(dec,artItEdge,articulationNode);
        }while(artItEdge != artFirstEdge);
    }

    MATRECfreeBlockArray(dec->env,&neighbourColors);
    return MATREC_OKAY;
}

static void intersectionOfAllPaths(MATRECGraphicDecomposition * dec, MATRECGraphicRowAddition *newRow,
                                   const reduced_member_id toCheck, int * const nodeNumPaths){
    int * intersectionPathDepth = newRow->intersectionPathDepth;
    spqr_node * intersectionPathParent = newRow->intersectionPathParent;

    //First do a dfs over the tree, storing all the tree-parents and depths for each node
    //TODO: maybe cache this tree and also update it so we can prevent this DFS call?


    //pick an arbitrary node as root; we just use the first cutEdge here
    {
        spqr_node root = findEdgeHead(dec, newRow->cutEdges[newRow->reducedMembers[toCheck].firstCutEdge].edge);
        DFSCallData *pathSearchCallStack = newRow->intersectionDFSData;

        assert(intersectionPathDepth[root] == -1);
        assert(intersectionPathParent[root] == SPQR_INVALID_NODE);

        int pathSearchCallStackSize = 0;

        intersectionPathDepth[root] = 0;
        intersectionPathParent[root] = SPQR_INVALID_NODE;

        pathSearchCallStack[0].node = root;
        pathSearchCallStack[0].nodeEdge = getFirstNodeEdge(dec, root);
        pathSearchCallStackSize++;
        while (pathSearchCallStackSize > 0) {
            assert(pathSearchCallStackSize <= newRow->memIntersectionDFSData);
            DFSCallData *dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
            //cannot be a tree edge which is its parent
            if (edgeIsTree(dec, dfsData->nodeEdge) &&
                (pathSearchCallStackSize <= 1 ||
                 dfsData->nodeEdge != pathSearchCallStack[pathSearchCallStackSize - 2].nodeEdge)) {
                spqr_node head = findEdgeHeadNoCompression(dec, dfsData->nodeEdge);
                spqr_node tail = findEdgeTailNoCompression(dec, dfsData->nodeEdge);
                spqr_node other = head == dfsData->node ? tail : head;
                assert(other != dfsData->node);

                //We go up a level: add new node to the call stack
                pathSearchCallStack[pathSearchCallStackSize].node = other;
                pathSearchCallStack[pathSearchCallStackSize].nodeEdge = getFirstNodeEdge(dec, other);
                //Every time a new node is discovered/added, we update its parent and depth information
                assert(intersectionPathDepth[other] == -1);
                assert(intersectionPathParent[other] == SPQR_INVALID_NODE);
                intersectionPathParent[other] = dfsData->node;
                intersectionPathDepth[other] = pathSearchCallStackSize;
                ++pathSearchCallStackSize;
                continue;
            }
            do {
                dfsData->nodeEdge = getNextNodeEdge(dec, dfsData->nodeEdge, dfsData->node);
                if (dfsData->nodeEdge == getFirstNodeEdge(dec, dfsData->node)) {
                    --pathSearchCallStackSize;
                    dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
                } else {
                    break;
                }
            } while (pathSearchCallStackSize > 0);
        }
    }

    //For each cut edge, trace back both ends until they meet
    cut_edge_id cutEdge = newRow->reducedMembers[toCheck].firstCutEdge;
    do{
        spqr_edge edge = newRow->cutEdges[cutEdge].edge;
        cutEdge = newRow->cutEdges[cutEdge].nextMember;

        //Iteratively jump up to the parents until they reach a common parent
        spqr_node source = findEdgeHead(dec, edge);
        spqr_node target = findEdgeTail(dec, edge);
        int sourceDepth = intersectionPathDepth[source];
        int targetDepth = intersectionPathDepth[target];
        nodeNumPaths[source]++;
        nodeNumPaths[target]++;

        while (sourceDepth > targetDepth){
            assert(source != target);
            source = intersectionPathParent[source];
            nodeNumPaths[source]++;
            --sourceDepth;
        }
        while(targetDepth > sourceDepth){
            assert(source != target);
            target = intersectionPathParent[target];
            nodeNumPaths[target]++;
            --targetDepth;
        }
        while(source != target && targetDepth >= 0){
            source = intersectionPathParent[source];
            target = intersectionPathParent[target];
            nodeNumPaths[source]++;
            nodeNumPaths[target]++;
            --targetDepth;
        }
        //In all the above, the lowest common ancestor is increased twice, so we correct for it ad-hoc
        nodeNumPaths[source]--;
        assert(SPQRnodeIsValid(source) && SPQRnodeIsValid(target));
        assert(source == target);

    }while (cutEdgeIsValid(cutEdge));

}

static void addArticulationNode(MATRECGraphicRowAddition *newRow, spqr_node articulationNode){
#ifndef NDEBUG
    for (int i = 0; i < newRow->numArticulationNodes; ++i) {
        assert(newRow->articulationNodes[i] != articulationNode);
    }
#endif
    newRow->articulationNodes[newRow->numArticulationNodes] = articulationNode;
    ++newRow->numArticulationNodes;
}
static void articulationPoints(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition * newRow, ArticulationNodeInformation *nodeInfo, reduced_member_id reducedMember){
    const bool * edgeRemoved = newRow->isEdgeCut;

    int rootChildren = 0;
    spqr_node root_node = findEdgeHead(dec, getFirstMemberEdge(dec, newRow->reducedMembers[reducedMember].member));;

    ArticulationPointCallStack * callStack = newRow->artDFSData;

    int depth = 0;
    int time = 1;

    callStack[depth].edge = getFirstNodeEdge(dec,root_node);
    callStack[depth].node = root_node;
    callStack[depth].parent = SPQR_INVALID_NODE;
    callStack[depth].isAP = false;

    nodeInfo[root_node].low = time;
    nodeInfo[root_node].discoveryTime = time;

    while(depth >= 0){
        if(!edgeRemoved[callStack[depth].edge]){
            spqr_node node = callStack[depth].node;
            spqr_node head = findEdgeHead(dec, callStack[depth].edge);
            spqr_node tail = findEdgeTail(dec, callStack[depth].edge);
            spqr_node otherNode = node == head ? tail : head;
            if(otherNode != callStack[depth].parent){
                if(nodeInfo[otherNode].discoveryTime == 0){
                    if(depth == 0){
                        rootChildren++;
                    }
                    //recursive call
                    ++depth;
                    assert(depth < newRow->memArtDFSData);
                    callStack[depth].parent = node;
                    callStack[depth].node = otherNode;
                    callStack[depth].edge = getFirstNodeEdge(dec,otherNode);
                    callStack[depth].isAP = false;

                    ++time;
                    nodeInfo[otherNode].low = time;
                    nodeInfo[otherNode].discoveryTime = time;
                    continue;

                }else{
                    nodeInfo[node].low = min(nodeInfo[node].low, nodeInfo[otherNode].discoveryTime);
                }
            }
        }

        while(true){
            callStack[depth].edge = getNextNodeEdge(dec,callStack[depth].edge,callStack[depth].node);
            if(callStack[depth].edge != getFirstNodeEdge(dec,callStack[depth].node)) break;
            --depth;
            if (depth < 0) break;

            spqr_node current_node = callStack[depth].node;
            spqr_node other_node = callStack[depth + 1].node;
            nodeInfo[current_node].low = min(nodeInfo[current_node].low,
                                             nodeInfo[other_node].low);
            if (depth != 0 &&
                !callStack[depth].isAP &&
                nodeInfo[current_node].discoveryTime <= nodeInfo[other_node].low) {
                addArticulationNode(newRow, current_node);
                callStack[depth].isAP = true;
            }
        }

    }
    if(rootChildren > 1 ){
        addArticulationNode(newRow,root_node);
    }
}

static void rigidConnectedColoringRecursive(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition * newRow, spqr_node articulationNode,
                                            spqr_node firstProcessNode, bool *isGood){
    const bool * isEdgeCut = newRow->isEdgeCut;
    COLOR_STATUS * nodeColors = newRow->nodeColors;
    ColorDFSCallData * data = newRow->colorDFSData;

    data[0].node = firstProcessNode;
    data[0].edge = getFirstNodeEdge(dec,firstProcessNode);
    newRow->nodeColors[firstProcessNode] = COLOR_FIRST;

    int depth = 0;
    while(depth >= 0){
        assert(depth < newRow->memColorDFSData);
        ColorDFSCallData * callData = &data[depth];
        spqr_node head = findEdgeHead(dec, callData->edge);
        spqr_node tail = findEdgeTail(dec, callData->edge);
        spqr_node otherNode = callData->node == head ? tail : head;
        if(otherNode != articulationNode){
            COLOR_STATUS currentColor = nodeColors[callData->node];
            COLOR_STATUS otherColor = nodeColors[otherNode];
            if(otherColor == UNCOLORED){
                if(isEdgeCut[callData->edge]){
                    nodeColors[otherNode] = currentColor == COLOR_FIRST ? COLOR_SECOND : COLOR_FIRST; //reverse the colors
                }else{
                    nodeColors[otherNode] = currentColor;
                }
                callData->edge = getNextNodeEdge(dec,callData->edge,callData->node);

                depth++;
                assert(depth < newRow->memColorDFSData);
                data[depth].node = otherNode;
                data[depth].edge = getFirstNodeEdge(dec,otherNode);
                continue;
            }else if(isEdgeCut[callData->edge] ^ (otherColor != currentColor)){
                *isGood = false;
                break;
            }
        }
        callData->edge = getNextNodeEdge(dec,callData->edge,callData->node);
        while(depth >= 0 && data[depth].edge == getFirstNodeEdge(dec,data[depth].node)){
            --depth;
        }
    }
}

static void rigidConnectedColoring(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                   const reduced_member_id reducedMember, const spqr_node node, bool * const isGood){

    spqr_node firstProcessNode;
    {
        spqr_edge edge = getFirstNodeEdge(dec, node);
        assert(SPQRedgeIsValid(edge));
        spqr_node head = findEdgeHead(dec, edge);
        spqr_node tail = findEdgeTail(dec, edge);
        //arbitrary way to select the first node
        firstProcessNode = head == node ? tail : head;

#ifndef NDEBUG
        {

            spqr_member member = findEdgeMemberNoCompression(dec, edge);
            spqr_edge firstEdge = getFirstMemberEdge(dec, member);
            spqr_edge memberEdge = firstEdge;
            do{
                assert(newRow->nodeColors[findEdgeHeadNoCompression(dec,memberEdge)] == UNCOLORED);
                assert(newRow->nodeColors[findEdgeTailNoCompression(dec,memberEdge)] == UNCOLORED);
                memberEdge = getNextMemberEdge(dec,memberEdge);
            }while(firstEdge != memberEdge);
        }
#endif
    }
    assert(SPQRnodeIsValid(firstProcessNode) && firstProcessNode != node);
    *isGood = true;
    newRow->reducedMembers[reducedMember].coloredNode = firstProcessNode;
    rigidConnectedColoringRecursive(dec,newRow,node,firstProcessNode,isGood);

    // Need to zero all colors for next attempts if we failed
    if(!(*isGood)){
        zeroOutColors(dec,newRow,firstProcessNode);
        newRow->reducedMembers[reducedMember].coloredNode = SPQR_INVALID_NODE;
    }else{
        zeroOutColorsExceptNeighbourhood(dec,newRow,node, firstProcessNode);
        newRow->reducedMembers[reducedMember].coloredNode = node;
    }
}

static spqr_node checkNeighbourColoringArticulationNode(MATRECGraphicDecomposition * dec, MATRECGraphicRowAddition *newRow,
                                                        const spqr_node articulationNode, spqr_edge * const adjacentSplittingEdge){
    spqr_node firstSideCandidate = SPQR_INVALID_NODE;
    spqr_node secondSideCandidate = SPQR_INVALID_NODE;
    spqr_edge firstSideEdge = SPQR_INVALID_EDGE;
    spqr_edge secondSideEdge = SPQR_INVALID_EDGE;
    int numFirstSide = 0;
    int numSecondSide = 0;

    spqr_edge firstEdge = getFirstNodeEdge(dec, articulationNode);
    spqr_edge moveEdge = firstEdge;
    do{
        spqr_node head = findEdgeHead(dec, moveEdge);
        spqr_node tail = findEdgeTail(dec, moveEdge);
        spqr_node otherNode = articulationNode == head ? tail : head;
        assert(newRow->nodeColors[otherNode] != UNCOLORED);
        //TODO: bit duplicate logic here? Maybe a nice way to fix?
        if((newRow->nodeColors[otherNode] == COLOR_FIRST) ^ newRow->isEdgeCut[moveEdge] ){
            if(numFirstSide == 0 && edgeIsTree(dec,moveEdge)){
                firstSideCandidate = otherNode;
                firstSideEdge = moveEdge;
            }else if (numFirstSide > 0){
                firstSideCandidate = SPQR_INVALID_NODE;
            }
            ++numFirstSide;
        }else{
            if(numSecondSide == 0 && edgeIsTree(dec,moveEdge)){
                secondSideCandidate = otherNode;
                secondSideEdge = moveEdge;
            }else if (numSecondSide > 0){
                secondSideCandidate = SPQR_INVALID_NODE;
            }
            ++numSecondSide;
        }
        moveEdge = getNextNodeEdge(dec,moveEdge,articulationNode);
    }while(moveEdge != firstEdge);

    if(numFirstSide == 1){
        *adjacentSplittingEdge = firstSideEdge;
        return firstSideCandidate;
    }else if (numSecondSide == 1){
        *adjacentSplittingEdge = secondSideEdge;
        return secondSideCandidate;
    }
    return SPQR_INVALID_NODE;
}


static void rigidGetSplittableArticulationPointsOnPath(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                                       const reduced_member_id toCheck, NodePair * const pair){
    int totalNumNodes = getNumNodes(dec);
    int * nodeNumPaths = newRow->crossingPathCount;

    for (int i = 0; i < totalNumNodes; ++i) {
        nodeNumPaths[i] = 0;
    }

    intersectionOfAllPaths(dec,newRow,toCheck,nodeNumPaths);

    newRow->numArticulationNodes = 0;

    ArticulationNodeInformation * artNodeInfo = newRow->articulationNodeSearchInfo;
    for (int i = 0; i < totalNumNodes; ++i) { //clean up can not easily be done in the search, unfortunately
        artNodeInfo[i].low = 0 ;
        artNodeInfo[i].discoveryTime = 0;
    }

    articulationPoints(dec,newRow,artNodeInfo,toCheck);

    int numCutEdges = newRow->reducedMembers[toCheck].numCutEdges;
    NodePairEmptyInitialize(&newRow->reducedMembers[toCheck].splitting_nodes);
    for (int i = 0; i < newRow->numArticulationNodes; i++) {
        spqr_node articulationNode = newRow->articulationNodes[i];
        assert(nodeIsRepresentative(dec, articulationNode));
        bool isOnPath = nodeNumPaths[articulationNode] == numCutEdges;
        if (isOnPath &&
            (NodePairIsEmpty(pair) || pair->first == articulationNode || pair->second == articulationNode)) {
            bool isGood = true;
            rigidConnectedColoring(dec, newRow, toCheck, articulationNode, &isGood);
            if (isGood) {
                NodePairInsert(&newRow->reducedMembers[toCheck].splitting_nodes, articulationNode);

                spqr_edge adjacentSplittingEdge = SPQR_INVALID_EDGE;
                spqr_node adjacentSplittingNode = checkNeighbourColoringArticulationNode(dec, newRow, articulationNode, &adjacentSplittingEdge);
                if (SPQRnodeIsValid(adjacentSplittingNode) &&
                    (NodePairIsEmpty(pair) || pair->first == adjacentSplittingNode ||
                     pair->second == adjacentSplittingNode)) {
                    bool isArticulationNode = false;
                    for (int j = 0; j < newRow->numArticulationNodes; ++j) {
                        if (newRow->articulationNodes[j] == adjacentSplittingNode) {
                            isArticulationNode = true;
                            break;
                        }
                    }
                    if (isArticulationNode) {
                        NodePairInsert(&newRow->reducedMembers[toCheck].splitting_nodes, adjacentSplittingNode);
                        assert(NodePairHasTwo(&newRow->reducedMembers[toCheck].splitting_nodes));
                        newRow->reducedMembers[toCheck].articulationEdge = adjacentSplittingEdge;
                        //Cleaning up the colors for next iterations...
                        {
                            spqr_edge firstNodeEdge = getFirstNodeEdge(dec, articulationNode);
                            spqr_edge itEdge = firstNodeEdge;
                            do {
                                spqr_node head = findEdgeHead(dec, itEdge);
                                spqr_node tail = findEdgeTail(dec, itEdge);
                                spqr_node otherNode = articulationNode == head ? tail : head;
                                newRow->nodeColors[otherNode] = UNCOLORED;
                                itEdge = getNextNodeEdge(dec, itEdge, articulationNode);
                            } while (itEdge != firstNodeEdge);

                        }
                    }
                }
                break;
            }
        }
    }

}

static void rowDetermineTypeRigid(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                               const reduced_member_id toCheckMember, const spqr_edge markerToOther,
                               const reduced_member_id otherMember, const spqr_edge markerToCheck){
    assert(newRow->reducedMembers[toCheckMember].numCutEdges > 0);//Checking for propagation only makes sense if there is at least one cut edge
    NodePair * nodePair = &newRow->reducedMembers[toCheckMember].splitting_nodes;
    {
        spqr_node head = findEdgeHead(dec, markerToOther);
        spqr_node tail = findEdgeTail(dec, markerToOther);
        NodePairInitialize(nodePair,head,tail);
    }
    if(newRow->reducedMembers[toCheckMember].numCutEdges == 1){
        //If there is a single cut edge, the two splitting points are its ends and thus there is no propagation

        spqr_edge cutEdge = newRow->cutEdges[newRow->reducedMembers[toCheckMember].firstCutEdge].edge;
        spqr_node head = findEdgeHead(dec, cutEdge);
        spqr_node tail = findEdgeTail(dec, cutEdge);
        NodePairIntersection(nodePair, head, tail);

        if(!NodePairIsEmpty(nodePair)){
            newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
            newRow->reducedMembers[toCheckMember].allHaveCommonNode = true;
        }else{
            newRow->reducedMembers[toCheckMember].type = TYPE_NOT_GRAPHIC;
            newRow->remainsGraphic = false;
        }
        assert(!NodePairHasTwo(nodePair)); //if this is not the case, there is a parallel edge
        return;
    }
    assert(newRow->reducedMembers[toCheckMember].numCutEdges > 1);
    spqr_edge articulationEdge = SPQR_INVALID_EDGE;
    NodePair pair = rigidDetermineAllAdjacentSplittableNodes(dec,newRow,toCheckMember,&articulationEdge);
    if(!NodePairIsEmpty(&pair)){
        if(SPQRedgeIsValid(articulationEdge)){
            if(articulationEdge == markerToOther){
                assert(edgeIsTree(dec,markerToOther));
                assert(NodePairHasTwo(&pair));
                newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
                createCutEdge(newRow,markerToCheck,otherMember);
                return;
            }
            NodePairIntersection(nodePair,pair.first,pair.second);
            assert(!NodePairHasTwo(nodePair)); // graph was not simple
            if(!NodePairIsEmpty(nodePair)){
                newRow->reducedMembers[toCheckMember].allHaveCommonNode = true;
                if(nodePair->first == pair.second){
                    newRow->reducedMembers[toCheckMember].articulationEdge = articulationEdge;
                }

                newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
            }else{
                newRow->reducedMembers[toCheckMember].type = TYPE_NOT_GRAPHIC;
                newRow->remainsGraphic = false;
            }

            return;
        }
        assert(!NodePairHasTwo(&pair));
        NodePairIntersection(nodePair,pair.first,pair.second);
        if(NodePairIsEmpty(nodePair)){
            newRow->reducedMembers[toCheckMember].type = TYPE_NOT_GRAPHIC;
            newRow->remainsGraphic = false;
        }else{
            newRow->reducedMembers[toCheckMember].allHaveCommonNode = true;
            newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
        }
        return;
    }

    NodePair copy = newRow->reducedMembers[toCheckMember].splitting_nodes; ; //At this point, this is simply the two nodes which were adjacent
    NodePairEmptyInitialize(nodePair);
    //TODO: below function mutates this internally... use a cleaner interface
    rigidGetSplittableArticulationPointsOnPath(dec,newRow,toCheckMember,&copy);

    if(NodePairIsEmpty(nodePair)){
        newRow->reducedMembers[toCheckMember].type = TYPE_NOT_GRAPHIC;
        newRow->remainsGraphic = false;
        return;
    }else if (NodePairHasTwo(nodePair)){
        assert(findEdgeHead(dec,markerToOther) == nodePair->first || findEdgeHead(dec,markerToOther) == nodePair->second);
        assert(findEdgeTail(dec,markerToOther) == nodePair->first || findEdgeTail(dec,markerToOther) == nodePair->second);
        newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
        createCutEdge(newRow,markerToCheck,otherMember);
        return;
    }else{
        newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
        return;
    }
}

static RowReducedMemberType determineType(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                       const reduced_member_id toCheckMember, const spqr_edge markerToOther,
                                       const reduced_member_id otherMember, const spqr_edge markerToCheck){
    assert(newRow->reducedMembers[toCheckMember].type == TYPE_UNDETERMINED);
    switch (getMemberType(dec,newRow->reducedMembers[toCheckMember].member)) {
        case SPQR_MEMBERTYPE_RIGID:
        {
            rowDetermineTypeRigid(dec,newRow,toCheckMember,markerToOther,otherMember,markerToCheck);
            break;

        }
        case SPQR_MEMBERTYPE_PARALLEL:
        {
            //we can only propagate if the marker edge is a tree edge and all other edges are cut
            if(edgeIsTree(dec,markerToOther) &&
               newRow->reducedMembers[toCheckMember].numCutEdges == (getNumMemberEdges(dec,newRow->reducedMembers[toCheckMember].member) - 1)){
                newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
                createCutEdge(newRow,markerToCheck,otherMember); //TODO: remove old cut edges? Or reuse memory maybe?
            }else{
                //In all other cases, the bond can be split so that the result will be okay!
                newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
            }
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:
        {
            //Propagation only calls this function if the edge is tree already, so we do not check it here.
            assert(edgeIsTree(dec,markerToOther));
            assert(newRow->reducedMembers[toCheckMember].numCutEdges == 1);
            newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
            createCutEdge(newRow,markerToCheck,otherMember); //TODO: remove old cut edges? Or reuse memory maybe?
            break;
        }
        default:
            assert(false);
            newRow->reducedMembers[toCheckMember].type = TYPE_NOT_GRAPHIC;
    }

    return newRow->reducedMembers[toCheckMember].type;
}

static void propagateComponents(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow){
    int leafArrayIndex = 0;

    reduced_member_id leaf;
    reduced_member_id next;

    while(leafArrayIndex != newRow->numLeafMembers){
        leaf = newRow->leafMembers[leafArrayIndex];
        //next is invalid if the member is not in the reduced decomposition.
        next = newRow->reducedMembers[leaf].parent;
        spqr_edge parentMarker = markerToParent(dec, newRow->reducedMembers[leaf].member);
        if(next != INVALID_REDUCED_MEMBER && edgeIsTree(dec,parentMarker)){
            assert(reducedMemberIsValid(next));
            assert(SPQRedgeIsValid(parentMarker));
            RowReducedMemberType type = determineType(dec,newRow,leaf,parentMarker,next,markerOfParent(dec,newRow->reducedMembers[leaf].member));
            if(type == TYPE_PROPAGATED){
                ++newRow->reducedMembers[next].numPropagatedChildren;
                if(newRow->reducedMembers[next].numPropagatedChildren == newRow->reducedMembers[next].numChildren){
                    newRow->leafMembers[leafArrayIndex] = next;
                }else{
                    ++leafArrayIndex;
                }
            }else if(type == TYPE_NOT_GRAPHIC){
                return;
            }else{
                assert(type == TYPE_MERGED);
                ++leafArrayIndex;
            }
        }else{
            ++leafArrayIndex;
        }

    }

    for (int j = 0; j < newRow->numReducedComponents; ++j) {
        //The reduced root might be a leaf as well: we propagate it last
        reduced_member_id root = newRow->reducedComponents[j].root;

        while(true){
            if(newRow->reducedMembers[root].numPropagatedChildren == newRow->reducedMembers[root].numChildren -1){
                //TODO: bit ugly, have to do a linear search for the child
                reduced_member_id child = INVALID_REDUCED_MEMBER;
                spqr_edge markerToChild = SPQR_INVALID_EDGE;
                for (children_idx i = newRow->reducedMembers[root].firstChild; i < newRow->reducedMembers[root].firstChild + newRow->reducedMembers[root].numChildren; ++i) {
                    reduced_member_id childReduced = newRow->childrenStorage[i];
                    if(newRow->reducedMembers[childReduced].type != TYPE_PROPAGATED){
                        child = childReduced;
                        markerToChild = markerOfParent(dec,newRow->reducedMembers[child].member);
                        break;
                    }
                }
                assert(SPQRmemberIsValid(newRow->reducedMembers[child].member));
                assert(SPQRedgeIsValid(markerToChild));
                if(!edgeIsTree(dec,markerToChild)){
                    break;
                }
                RowReducedMemberType type = determineType(dec,newRow,root,markerToChild,child,markerToParent(dec,newRow->reducedMembers[child].member));
                if(type == TYPE_PROPAGATED){
                    root = child;
                }else if(type == TYPE_NOT_GRAPHIC){
                    return;
                }else{
                    break;
                }
            }else{
                break;
            }
        }
        newRow->reducedComponents[j].root = root;
        newRow->reducedMembers[root].parent = INVALID_REDUCED_MEMBER;
    }
}

static NodePair
rigidDetermineCandidateNodesFromAdjacentComponents(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow, const reduced_member_id toCheck) {

    NodePair pair;
    NodePairEmptyInitialize(&pair);
    NodePair * nodePair = &pair;
    assert(NodePairIsEmpty(nodePair));

    //take union of children's edges nodes to find one or two candidates
    for (children_idx i = newRow->reducedMembers[toCheck].firstChild;
         i < newRow->reducedMembers[toCheck].firstChild + newRow->reducedMembers[toCheck].numChildren; ++i) {
        reduced_member_id reducedChild = newRow->childrenStorage[i];
        if (newRow->reducedMembers[reducedChild].type != TYPE_PROPAGATED) {
            spqr_edge edge = markerOfParent(dec, newRow->reducedMembers[reducedChild].member);
            spqr_node head = findEdgeHead(dec, edge);
            spqr_node tail = findEdgeTail(dec, edge);
            if(NodePairIsEmpty(nodePair)){
                NodePairInitialize(nodePair,head,tail);
            }else{
                NodePairIntersection(nodePair, head, tail);
            }
            if (NodePairIsEmpty(nodePair)) {
                return pair;
            }
        }
    }
    if (reducedMemberIsValid(newRow->reducedMembers[toCheck].parent) &&
        newRow->reducedMembers[newRow->reducedMembers[toCheck].parent].type != TYPE_PROPAGATED) {

        spqr_edge edge = markerToParent(dec, newRow->reducedMembers[toCheck].member);
        spqr_node head = findEdgeHead(dec, edge);
        spqr_node tail = findEdgeTail(dec, edge);
        if(NodePairIsEmpty(nodePair)){
            NodePairInitialize(nodePair,head,tail);
        }else{
            NodePairIntersection(nodePair, head, tail);
        }
    }
    return pair;
}

static void determineTypeRigidMerging(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow, const reduced_member_id toCheck){
    NodePair * nodePair = &newRow->reducedMembers[toCheck].splitting_nodes;
    bool hasNoAdjacentMarkers = (newRow->reducedMembers[toCheck].numChildren - newRow->reducedMembers[toCheck].numPropagatedChildren) == 0 &&
                                (reducedMemberIsInvalid(newRow->reducedMembers[toCheck].parent) ||
                                 newRow->reducedMembers[newRow->reducedMembers[toCheck].parent].type == TYPE_PROPAGATED);
    //if the component is free standing
    if(hasNoAdjacentMarkers){
        assert(newRow->reducedMembers[toCheck].numCutEdges> 0);
        if(newRow->reducedMembers[toCheck].numCutEdges == 1){
            //TODO: with two splitting nodes maybe store the cut edge?
            newRow->reducedMembers[toCheck].type = TYPE_MERGED;
        }else{

            //take union of edge ends
            cut_edge_id cutEdgeIdx = newRow->reducedMembers[toCheck].firstCutEdge;
            spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
            spqr_node head = findEdgeHead(dec, cutEdge);
            spqr_node tail = findEdgeTail(dec, cutEdge);
            NodePairInitialize(nodePair,head,tail);

            while(cutEdgeIsValid(newRow->cutEdges[cutEdgeIdx].nextMember)){
                cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextMember;
                cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
                head = findEdgeHead(dec, cutEdge);
                tail = findEdgeTail(dec, cutEdge);
                NodePairIntersection(nodePair, head, tail);

                if(NodePairIsEmpty(nodePair)){
                    break;
                }
            }
            //Since there is two ore more edges, there can be at most one splitting node:
            assert(!NodePairHasTwo(nodePair));
            if(!NodePairIsEmpty(nodePair)){
                //All cut edges are adjacent to one node; either this node can be 'extended' or if the number of cut edges == degree of node -1,
                //the other edge is placed in series with the new row edge
                newRow->reducedMembers[toCheck].allHaveCommonNode = true;
                newRow->reducedMembers[toCheck].type = TYPE_MERGED;
                return;
            }
            NodePair emptyPair;
            NodePairEmptyInitialize(&emptyPair);
            rigidGetSplittableArticulationPointsOnPath(dec,newRow,toCheck,&emptyPair);
            if(NodePairIsEmpty(nodePair)){
                newRow->reducedMembers[toCheck].type = TYPE_NOT_GRAPHIC;
                newRow->remainsGraphic = false;
                return;
            }else{
                newRow->reducedMembers[toCheck].type = TYPE_MERGED;
                return;
            }
        }
        return;
    }

    NodePair adjacentMarkers = rigidDetermineCandidateNodesFromAdjacentComponents(dec,newRow,toCheck);

    if(NodePairIsEmpty(&adjacentMarkers)){
        NodePairEmptyInitialize(&newRow->reducedMembers[toCheck].splitting_nodes);
        newRow->reducedMembers[toCheck].type = TYPE_NOT_GRAPHIC;
        newRow->remainsGraphic = false;
        return;
    }
    //Check if all edges are adjacent or have an articulation point
    if(newRow->reducedMembers[toCheck].numCutEdges > 0){
        spqr_edge articulationEdgeToSecond = SPQR_INVALID_EDGE;
        NodePair edgeAdjacentPair = rigidDetermineAllAdjacentSplittableNodes(dec,newRow,toCheck,&articulationEdgeToSecond);
        if(!NodePairIsEmpty(&edgeAdjacentPair)){
            NodePairIntersection(&adjacentMarkers,edgeAdjacentPair.first,edgeAdjacentPair.second);
            if(NodePairIsEmpty(&adjacentMarkers)){
                NodePairEmptyInitialize(&newRow->reducedMembers[toCheck].splitting_nodes);
                newRow->reducedMembers[toCheck].type = TYPE_NOT_GRAPHIC;
                newRow->remainsGraphic = false;
            }else{
                newRow->reducedMembers[toCheck].allHaveCommonNode = true;
                assert(!NodePairHasTwo(&adjacentMarkers)); //graph should have been simple
                if(SPQRedgeIsValid(articulationEdgeToSecond) && adjacentMarkers.first == edgeAdjacentPair.second){
                    newRow->reducedMembers[toCheck].articulationEdge = articulationEdgeToSecond;
                }
                newRow->reducedMembers[toCheck].type = TYPE_MERGED;
                NodePairEmptyInitialize(&newRow->reducedMembers[toCheck].splitting_nodes);
                NodePairInsert(&newRow->reducedMembers[toCheck].splitting_nodes,adjacentMarkers.first);
            }
            return;
        }


        assert(NodePairIsEmpty(nodePair));
        rigidGetSplittableArticulationPointsOnPath(dec,newRow,toCheck,&adjacentMarkers);
        if(NodePairIsEmpty(nodePair)){
            newRow->reducedMembers[toCheck].type = TYPE_NOT_GRAPHIC;
            newRow->remainsGraphic = false;
            return;
        }
        assert(!NodePairHasTwo(nodePair)); //Graph was not simple
        newRow->reducedMembers[toCheck].type = TYPE_MERGED;
        return;
    }

    //No cut edges: simply take the point of the adjacent markers
    assert(!NodePairHasTwo(&adjacentMarkers));
    NodePairEmptyInitialize(&newRow->reducedMembers[toCheck].splitting_nodes);
    NodePairInsert(&newRow->reducedMembers[toCheck].splitting_nodes,adjacentMarkers.first);
    newRow->reducedMembers[toCheck].type = TYPE_MERGED;
}

static RowReducedMemberType determineTypeMerging(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow, const reduced_member_id toCheck){
    switch(getMemberType(dec,newRow->reducedMembers[toCheck].member)){
        case SPQR_MEMBERTYPE_RIGID:
        {
            determineTypeRigidMerging(dec,newRow,toCheck);
            break;
        }
        case SPQR_MEMBERTYPE_LOOP:
        case SPQR_MEMBERTYPE_PARALLEL:
        {
            newRow->reducedMembers[toCheck].type = TYPE_MERGED;
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:
        {
            int numNonPropagatedAdjacent = newRow->reducedMembers[toCheck].numChildren-newRow->reducedMembers[toCheck].numPropagatedChildren;
            if(reducedMemberIsValid(newRow->reducedMembers[toCheck].parent) &&
               newRow->reducedMembers[newRow->reducedMembers[toCheck].parent].type != TYPE_PROPAGATED){
                ++numNonPropagatedAdjacent;
            }

            assert(numNonPropagatedAdjacent != 1);
            if(numNonPropagatedAdjacent > 2){
                newRow->reducedMembers[toCheck].type = TYPE_NOT_GRAPHIC;
                newRow->remainsGraphic = false;
            }else{
                newRow->reducedMembers[toCheck].type = TYPE_MERGED;
            }
            break;
        }
        default:
        {
            assert(false);
            newRow->reducedMembers[toCheck].type = TYPE_NOT_GRAPHIC;
            break;
        }
    }
    return newRow->reducedMembers[toCheck].type;
}

static MATREC_ERROR allocateTreeSearchMemory(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow){
    int necessarySpace = newRow->numReducedMembers;
    if( necessarySpace > newRow->memMergeTreeCallData ){
        newRow->memMergeTreeCallData = max(2*newRow->memMergeTreeCallData,necessarySpace);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->mergeTreeCallData,(size_t) newRow->memMergeTreeCallData));
    }
    return MATREC_OKAY;
}

static void determineMergeableTypes(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow, reduced_member_id root){
    assert(newRow->numReducedMembers <= newRow->memMergeTreeCallData);

    int depth = 0;
    MergeTreeCallData * stack = newRow->mergeTreeCallData;

    stack[0].currentChild = newRow->reducedMembers[root].firstChild;
    stack[0].id = root;

    //First determine type of all children, then determine type of the node itself
    do{
        if(stack[depth].currentChild == newRow->reducedMembers[stack[depth].id].firstChild +
                                        newRow->reducedMembers[stack[depth].id].numChildren){
            if(newRow->reducedMembers[stack[depth].id].type == TYPE_UNDETERMINED){
                RowReducedMemberType type = determineTypeMerging(dec,newRow,stack[depth].id);
                if(type == TYPE_NOT_GRAPHIC){
                    assert(!newRow->remainsGraphic);
                    break;
                }

            }
            --depth;
            continue; //does this break when necessary?
        }
        reduced_member_id reducedChild = newRow->childrenStorage[stack[depth].currentChild];
        ++stack[depth].currentChild;
        if (newRow->reducedMembers[reducedChild].type != TYPE_PROPAGATED) {
            ++depth;
            assert(depth < newRow->memMergeTreeCallData);
            stack[depth].id = reducedChild;
            stack[depth].currentChild = newRow->reducedMembers[reducedChild].firstChild;
        }
    }while(depth >= 0);
}

static void cleanUpRowMemberInformation(MATRECGraphicRowAddition * newRow){
    //This loop is at the end as memberInformation is also used to assign the cut edges during propagation
    //Clean up the memberInformation array
    for (int i = 0; i < newRow->numReducedMembers; ++i) {
        newRow->memberInformation[newRow->reducedMembers[i].member].reducedMember = INVALID_REDUCED_MEMBER;
    }
#ifndef NDEBUG
    for (int i = 0; i < newRow->memMemberInformation; ++i) {
        assert(reducedMemberIsInvalid(newRow->memberInformation[i].reducedMember));
    }
#endif
}

static MATREC_ERROR rigidTransformEdgeIntoCycle(MATRECGraphicDecomposition *dec,
                                              const spqr_member member,
                                              const spqr_edge edge,
                                              NewRowInformation * const newRowInformation){
    //If a cycle already exists, just expand it with the new edge.
    spqr_member markerCycleMember = SPQR_INVALID_MEMBER;
    if (edge == markerToParent(dec, member)) {
        spqr_member parentMember = findMemberParent(dec, member);
        if (getMemberType(dec, parentMember) == SPQR_MEMBERTYPE_SERIES) {
            markerCycleMember = parentMember;
        }
    } else if (edgeIsMarker(dec, edge)) {
        spqr_member childMember = findEdgeChildMember(dec, edge);
        if (getMemberType(dec, childMember) == SPQR_MEMBERTYPE_SERIES) {
            markerCycleMember = childMember;
        }
    }
    if (markerCycleMember != SPQR_INVALID_MEMBER) {
        newRowInformation->member = markerCycleMember;
        return MATREC_OKAY;
    }

    //create a new cycle
    spqr_member new_cycle;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &new_cycle));
    bool isTreeEdge = edgeIsTree(dec,edge);
    bool parentMoved = markerToParent(dec,member) == edge;
    moveEdgeToNewMember(dec,edge,member,new_cycle);
    spqr_node cutHead = findEdgeHead(dec, edge);
    spqr_node cutTail = findEdgeTail(dec, edge);
    clearEdgeHeadAndTail(dec,edge);

    spqr_edge ignore = SPQR_INVALID_EDGE;
    spqr_edge rigidMarker = SPQR_INVALID_EDGE;

    if(parentMoved){
        MATREC_CALL(createMarkerPairWithReferences(dec,new_cycle,member,!isTreeEdge,&ignore,&rigidMarker));
    }else{
        MATREC_CALL(createMarkerPairWithReferences(dec,member,new_cycle,isTreeEdge,&rigidMarker,&ignore));

    }
    setEdgeHeadAndTail(dec,rigidMarker,cutHead,cutTail);

    newRowInformation->member = new_cycle;

    return MATREC_OKAY;
}

static MATREC_ERROR transformSingleRigid(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                       const reduced_member_id reducedMember,
                                       const spqr_member member,
                                       NewRowInformation * const newRowInformation){
    if(newRow->reducedMembers[reducedMember].numCutEdges == 1){
        //Cut edge is propagated into a cycle with new edge
        spqr_edge cutEdge = newRow->cutEdges[newRow->reducedMembers[reducedMember].firstCutEdge].edge;
        MATREC_CALL(rigidTransformEdgeIntoCycle(dec,member,cutEdge,newRowInformation));
        return MATREC_OKAY;
    }
    assert(newRow->reducedMembers[reducedMember].numCutEdges > 1);
    assert(!NodePairIsEmpty(&newRow->reducedMembers[reducedMember].splitting_nodes));

    spqr_node splitNode = newRow->reducedMembers[reducedMember].splitting_nodes.first;

    //For now this code assumes that all nodes are adjacent to the split node!!
    if(newRow->reducedMembers[reducedMember].allHaveCommonNode){
        if (newRow->reducedMembers[reducedMember].numCutEdges != nodeDegree(dec, splitNode) - 1) {
            //Create a new node; move all cut edges end of split node to it and add new edge between new node and split node
            spqr_node newNode = SPQR_INVALID_NODE;
            MATREC_CALL(createNode(dec, &newNode));
            newRowInformation->member = member;
            newRowInformation->firstNode = newNode;
            newRowInformation->secondNode = splitNode;

            cut_edge_id cutEdgeIdx = newRow->reducedMembers[reducedMember].firstCutEdge;
            do {
                spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
                spqr_node edgeHead = findEdgeHead(dec, cutEdge);
                if (edgeHead == splitNode) {
                    changeEdgeHead(dec, cutEdge, edgeHead, newNode);
                } else {
                    changeEdgeTail(dec, cutEdge, findEdgeTail(dec, cutEdge), newNode);
                }

                cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextMember;
            } while (cutEdgeIsValid(cutEdgeIdx));

            return MATREC_OKAY;
        } else {
            //We need to find the non-cut edge. By definition, this must be a tree edge
            //TODO: searching for the edge is a bit ugly, but we need to do it somewhere

            spqr_edge firstEdge = getFirstNodeEdge(dec, splitNode);
            spqr_edge moveEdge = firstEdge;
            do {
                if (edgeIsTree(dec, moveEdge)) {
                    break;
                }
                moveEdge = getNextNodeEdge(dec, moveEdge, splitNode);
            } while (moveEdge != firstEdge);
            assert(edgeIsTree(dec, moveEdge));

            MATREC_CALL(rigidTransformEdgeIntoCycle(dec,member,moveEdge,newRowInformation));
            return MATREC_OKAY;
        }
    }


    if(NodePairHasTwo(&newRow->reducedMembers[reducedMember].splitting_nodes)){
        //Two splitting nodes: identify the edge and either make it a cycle or prolong it, if it is already a marker to a cycle
        spqr_edge makeCycleEdge = newRow->reducedMembers[reducedMember].articulationEdge;

        assert(edgeIsTree(dec,makeCycleEdge));
        assert(SPQRedgeIsValid(makeCycleEdge));

        MATREC_CALL(rigidTransformEdgeIntoCycle(dec,member,makeCycleEdge,newRowInformation));
        return MATREC_OKAY;
    }
    assert(!NodePairIsEmpty(&newRow->reducedMembers[reducedMember].splitting_nodes));
    //Multiple edges without a common end: need to use coloring information
    int numFirstColor = 0;
    int numSecondColor = 0;

    spqr_edge firstNodeEdge = getFirstNodeEdge(dec, splitNode);
    spqr_edge iterEdge = firstNodeEdge;
    do{
        spqr_node head = findEdgeHead(dec, iterEdge);
        spqr_node tail = findEdgeTail(dec, iterEdge);
        spqr_node other = head == splitNode ? tail : head;
        if(newRow->nodeColors[other] == COLOR_FIRST){
            numFirstColor++;
        }else{
            numSecondColor++;
        }
        iterEdge = getNextNodeEdge(dec,iterEdge,splitNode);
    }while(iterEdge != firstNodeEdge);

    COLOR_STATUS toNewNodeColor = numFirstColor < numSecondColor ? COLOR_FIRST : COLOR_SECOND; //TODO: sanity check this?

    spqr_node newNode = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec, &newNode));


    //TODO: fix duplication? (with where?)
    {
        firstNodeEdge = getFirstNodeEdge(dec,splitNode);
        iterEdge = firstNodeEdge;
        do{
            bool isCut = newRow->isEdgeCut[iterEdge];
            spqr_node otherHead = findEdgeHead(dec, iterEdge);
            spqr_node otherTail = findEdgeTail(dec, iterEdge);
            spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
            bool isMoveColor = newRow->nodeColors[otherEnd] == toNewNodeColor;
            spqr_edge nextEdge = getNextNodeEdge(dec, iterEdge, splitNode); //Need to do this before we modify the edge :)

            bool changeEdgeEnd = (isCut && isMoveColor) || (!isCut && !isMoveColor);
            if(changeEdgeEnd){
                if(otherHead == splitNode){
                    changeEdgeHead(dec,iterEdge,otherHead,newNode);
                }else{
                    changeEdgeTail(dec,iterEdge,otherTail,newNode);
                }
            }
            newRow->nodeColors[otherEnd] = UNCOLORED; //Clean up

            //Ugly hack to make sure we can iterate neighbourhood whilst changing edge ends.
            spqr_edge previousEdge = iterEdge;
            iterEdge = nextEdge;
            if(iterEdge == firstNodeEdge){
                break;
            }
            if(changeEdgeEnd && previousEdge == firstNodeEdge){
                firstNodeEdge = iterEdge;
            }
        }while(true);
        newRow->reducedMembers[reducedMember].coloredNode = SPQR_INVALID_NODE;
    }
    newRowInformation->member = member;
    newRowInformation->firstNode = splitNode;
    newRowInformation->secondNode = newNode;
    return MATREC_OKAY;
}


static MATREC_ERROR splitParallelRowAddition(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                const reduced_member_id reducedMember,
                                const spqr_member member,
                                spqr_member * const loopMember){
    assert(newRow->reducedMembers[reducedMember].numCutEdges > 0);

    int numCutEdges = newRow->reducedMembers[reducedMember].numCutEdges;
    int numParallelEdges = getNumMemberEdges(dec,member);

    bool createCutParallel = numCutEdges > 1;
    bool convertOriginalParallel = (numCutEdges + 1) == numParallelEdges;

    if(!createCutParallel && convertOriginalParallel){
        assert(getNumMemberEdges(dec,member) == 2);
        assert(newRow->reducedMembers[reducedMember].numCutEdges == 1);
        *loopMember = member;
        return MATREC_OKAY;
    }

    spqr_member cutMember = SPQR_INVALID_MEMBER;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

    cut_edge_id cutEdgeIdx = newRow->reducedMembers[reducedMember].firstCutEdge;
    assert(cutEdgeIsValid(cutEdgeIdx));
    bool parentCut = false;

    while(cutEdgeIsValid(cutEdgeIdx)){
        spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
        cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextMember;
        moveEdgeToNewMember(dec,cutEdge,member,cutMember);
        if (cutEdge == markerToParent(dec,member)){
            parentCut = true;
        }
    }
    if(convertOriginalParallel == createCutParallel){
        if(parentCut){
            MATREC_CALL(createMarkerPair(dec,cutMember,member,true));
        }else{
            MATREC_CALL(createMarkerPair(dec,member,cutMember,false));
        }
        *loopMember = convertOriginalParallel ? member : cutMember;
        return MATREC_OKAY;
    }

    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, loopMember));

    //TODO: check tree directions here...
    //TODO Can probably eliminate some branches here
    if(parentCut){
        MATREC_CALL(createMarkerPair(dec,cutMember,*loopMember,true));
        MATREC_CALL(createMarkerPair(dec,*loopMember,member,true));
    }else{
        MATREC_CALL(createMarkerPair(dec,member,*loopMember,false));
        MATREC_CALL(createMarkerPair(dec,*loopMember,cutMember,false));
    }
    return MATREC_OKAY;
}

static MATREC_ERROR transformSingleParallel(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                          const reduced_member_id reducedMember,
                                          const spqr_member member,
                                          spqr_member * const newRowMember){
    spqr_member loopMember = SPQR_INVALID_MEMBER;
    MATREC_CALL(splitParallelRowAddition(dec,newRow,reducedMember,member,&loopMember));

    //TODO: fix with adjacent cycle
    //Make loop member a series and elongate it with the new row,
    //The only exception is if the non-cut edge is a child or parent marker to a series. In that case, the series member is elongated and the loop is merged
    //TODO: now we iterate to find the non-cut edge, with better memory layout this would be much simpler

    bool adjacent_series_marker = false;
    {
        spqr_edge first_edge = getFirstMemberEdge(dec, loopMember);
        spqr_edge edge = first_edge;

        do {
            if(edgeIsMarker(dec,edge)){
                spqr_member child = findEdgeChildMember(dec, edge);
                if(getMemberType(dec,child) == SPQR_MEMBERTYPE_SERIES){
                    adjacent_series_marker = true;
                    *newRowMember = child;
                    break;
                }
            }else if(markerToParent(dec,member) == edge) {
                spqr_member parent = findMemberParent(dec, loopMember);
                if (getMemberType(dec, parent) == SPQR_MEMBERTYPE_SERIES) {
                    adjacent_series_marker = true;
                    *newRowMember = parent;
                    break;
                }
            }
            edge = getNextMemberEdge(dec, edge);
        } while (edge != first_edge);
    }
    if(!adjacent_series_marker){
        changeLoopToSeries(dec,loopMember);
        *newRowMember = loopMember;
    }else{
        mergeLoop(dec, loopMember);
    }



    return MATREC_OKAY;
}

/**
 * Splits a parallel member into two parallel members connected by a loop, based on which edges are cut.
 * For both of the bonds if they would have only 2 edges, they are merged into the middle bond
 * @return
 */
static MATREC_ERROR splitParallelMerging(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                       const reduced_member_id reducedMember,
                                       const spqr_member member,
                                       spqr_member * const pMergeMember,
                                       spqr_edge * const cutRepresentative){
    //When merging, we cannot have propagated members;
    assert(newRow->reducedMembers[reducedMember].numCutEdges < (getNumMemberEdges(dec,member)-1));

    int numMergeableAdjacent = newRow->reducedMembers[reducedMember].numChildren - newRow->reducedMembers[reducedMember].numPropagatedChildren;
    if(reducedMemberIsValid(newRow->reducedMembers[reducedMember].parent) &&
       newRow->reducedMembers[newRow->reducedMembers[reducedMember].parent].type == TYPE_MERGED){
        numMergeableAdjacent++;
    }

    int numCutEdges = newRow->reducedMembers[reducedMember].numCutEdges;
    //All edges which are not in the mergeable decomposition or cut
    int numBaseSplitAwayEdges = getNumMemberEdges(dec,member) - numMergeableAdjacent - numCutEdges ;

    bool createCutParallel = numCutEdges > 1;
    bool keepOriginalParallel = numBaseSplitAwayEdges  <= 1;

    bool parentCut = false;

    spqr_member cutMember = SPQR_INVALID_MEMBER;
    if(createCutParallel){
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

        cut_edge_id cutEdgeIdx = newRow->reducedMembers[reducedMember].firstCutEdge;
        assert(cutEdgeIsValid(cutEdgeIdx));

        while(cutEdgeIsValid(cutEdgeIdx)){
            spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
            cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextMember;
            moveEdgeToNewMember(dec,cutEdge,member,cutMember);
            if (cutEdge == markerToParent(dec,member)){
                parentCut = true;
            }
        }

    }else if(numCutEdges == 1){
        *cutRepresentative = newRow->cutEdges[newRow->reducedMembers[reducedMember].firstCutEdge].edge;
    }

    spqr_edge noCutRepresentative = SPQR_INVALID_EDGE;
    spqr_member mergingMember = member;
    bool parentToMergingMember = false;
    bool treeToMergingMember = false;
    if(!keepOriginalParallel){
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &mergingMember));
        //move all mergeable children and parent edges to the mergingMember
        for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
             i < newRow->reducedMembers[reducedMember].firstChild + newRow->reducedMembers[reducedMember].numChildren; ++i) {
            reduced_member_id child = newRow->childrenStorage[i];
            if(newRow->reducedMembers[child].type == TYPE_MERGED){
                spqr_edge moveEdge = markerOfParent(dec, newRow->reducedMembers[child].member);
                moveEdgeToNewMember(dec,moveEdge,member,mergingMember);
                if(edgeIsTree(dec,moveEdge)){
                    treeToMergingMember = true;
                }
            }
        }
        reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
        if(reducedMemberIsValid(parent) &&
           newRow->reducedMembers[parent].type == TYPE_MERGED && !parentCut){
            spqr_edge moveEdge = markerToParent(dec, member);
            moveEdgeToNewMember(dec,moveEdge,member,mergingMember);
            parentToMergingMember = true;
            if(edgeIsTree(dec,moveEdge)){
                treeToMergingMember = true;
            }
        }
        //If there is only one cut edge, we also move it.
        if(SPQRedgeIsValid(*cutRepresentative)){
            if(*cutRepresentative == markerToParent(dec,member)){
                parentToMergingMember = true;
            }
            moveEdgeToNewMember(dec,*cutRepresentative,member,mergingMember);
        }
    }
    //TODO: can probably reduce branching a bit here.
    if(createCutParallel){
        spqr_edge ignoreArgument = SPQR_INVALID_EDGE;
        if(parentCut){
            MATREC_CALL(createMarkerPairWithReferences(dec,cutMember,mergingMember,true,&ignoreArgument,cutRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,mergingMember,cutMember,false,cutRepresentative,&ignoreArgument));
        }
    }
    if(!keepOriginalParallel){
        spqr_edge ignoreArgument = SPQR_INVALID_EDGE;
        if(parentToMergingMember){
            MATREC_CALL(createMarkerPairWithReferences(dec,mergingMember,member,!treeToMergingMember,&ignoreArgument,&noCutRepresentative));
        }else{
            if(parentCut){
                MATREC_CALL(createMarkerPairWithReferences(dec,mergingMember,member,!treeToMergingMember,&ignoreArgument,&noCutRepresentative));
            }else{
                MATREC_CALL(createMarkerPairWithReferences(dec,member,mergingMember,treeToMergingMember,&noCutRepresentative,&ignoreArgument));
            }
        }
    }
    *pMergeMember = mergingMember;
    return MATREC_OKAY;
}

static MATREC_ERROR splitSeriesMergingRowAddition(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                     const reduced_member_id reducedMember,
                                     const spqr_member member,
                                     spqr_member * const mergingMember, bool * const isCut){
    assert(getNumMemberEdges(dec,member) >= 3);
    * isCut = newRow->reducedMembers[reducedMember].numCutEdges > 0;

    if(getNumMemberEdges(dec,member) == 3){

        *mergingMember = member;
        return MATREC_OKAY;
    }
    //Split off the relevant part of the series member
    spqr_member mergingSeries = SPQR_INVALID_MEMBER;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &mergingSeries));
    //Move all marker edges which point to another component in the reduced decomposition to the new member
    //This should be exactly 2, as with 3 the result is not graphic anymore
    //move all mergeable children and parent edges to the mergingMember

    bool coTreeToMergingMember = false;
    bool parentToMergingMember = false;
    for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
         i < newRow->reducedMembers[reducedMember].firstChild + newRow->reducedMembers[reducedMember].numChildren; ++i) {
        reduced_member_id child = newRow->childrenStorage[i];
        if(newRow->reducedMembers[child].type == TYPE_MERGED){
            spqr_edge moveEdge = markerOfParent(dec, newRow->reducedMembers[child].member);
            moveEdgeToNewMember(dec,moveEdge,member,mergingSeries);
            if(!edgeIsTree(dec,moveEdge)){
                coTreeToMergingMember = true;
            }
        }
    }

    reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
    if(reducedMemberIsValid(parent) &&
       newRow->reducedMembers[parent].type == TYPE_MERGED ){
        spqr_edge moveEdge = markerToParent(dec, member);
        moveEdgeToNewMember(dec,moveEdge,member,mergingSeries);
        parentToMergingMember = true;
        if(!edgeIsTree(dec,moveEdge)){
            coTreeToMergingMember = true;
        }
    }
    if(parentToMergingMember){
        MATREC_CALL(createMarkerPair(dec,mergingSeries,member,coTreeToMergingMember));
    }else{
        MATREC_CALL(createMarkerPair(dec,member,mergingSeries,!coTreeToMergingMember));
    }

    *mergingMember = mergingSeries;
    assert(getNumMemberEdges(dec,mergingSeries) == 3 );
    assert(getNumMemberEdges(dec,member) >= 3);
    return MATREC_OKAY;
}

static MATREC_ERROR transformRoot(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                const reduced_member_id reducedMember,
                                NewRowInformation * const newRowInfo){
    assert(newRow->reducedMembers[reducedMember].type == TYPE_MERGED);

    spqr_member member = newRow->reducedMembers[reducedMember].member;

    //For any series or parallel member, 'irrelevant parts' are first split off into separate parallel and series members
    //For series and parallel members, we need to create nodes and correctly split them

    //TODO: split up into functions
    switch(getMemberType(dec,member)){
        case SPQR_MEMBERTYPE_RIGID: {
            spqr_node newNode = SPQR_INVALID_NODE;
            MATREC_CALL(createNode(dec,&newNode));
            //we should have identified a unique splitting node
            assert(!NodePairIsEmpty(&newRow->reducedMembers[reducedMember].splitting_nodes)
                   && !NodePairHasTwo(&newRow->reducedMembers[reducedMember].splitting_nodes));
            spqr_node splitNode = newRow->reducedMembers[reducedMember].splitting_nodes.first;

            newRowInfo->firstNode = splitNode;
            newRowInfo->secondNode = newNode;
            newRowInfo->member = member;

            if(newRow->reducedMembers[reducedMember].numCutEdges == 0){
                break;
            }
            if (newRow->reducedMembers[reducedMember].allHaveCommonNode) {
                if(SPQRedgeIsValid(newRow->reducedMembers[reducedMember].articulationEdge)){
                    spqr_edge articulationEdge =newRow->reducedMembers[reducedMember].articulationEdge;
                    spqr_node cutHead = findEdgeHead(dec, articulationEdge);
                    spqr_node cutTail = findEdgeTail(dec, articulationEdge);
                    if(cutHead == splitNode){
                        changeEdgeHead(dec,articulationEdge,cutHead,newNode);
                    }else{
                        assert(cutTail == splitNode);
                        changeEdgeTail(dec,articulationEdge,cutTail,newNode);
                    }
                }else{

                    cut_edge_id cutEdgeIdx = newRow->reducedMembers[reducedMember].firstCutEdge;
                    do{
                        spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
                        spqr_node cutHead = findEdgeHead(dec, cutEdge);
                        spqr_node cutTail = findEdgeTail(dec, cutEdge);
                        bool moveHead = cutHead == splitNode;
                        assert((moveHead && cutHead == splitNode) || (!moveHead && cutTail == splitNode));
                        if(moveHead){
                            changeEdgeHead(dec,cutEdge,cutHead,newNode);
                        }else{
                            changeEdgeTail(dec,cutEdge,cutTail,newNode);
                        }
                        cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextMember;
                    }while(cutEdgeIsValid(cutEdgeIdx));
                }

            } else {
                //TODO: fix duplication here.
                int numFirstColor = 0;
                int numSecondColor = 0;

                spqr_edge firstNodeEdge = getFirstNodeEdge(dec, splitNode);
                spqr_edge iterEdge = firstNodeEdge;
                do{
                    spqr_node head = findEdgeHead(dec, iterEdge);
                    spqr_node tail = findEdgeTail(dec, iterEdge);
                    spqr_node other = head == splitNode ? tail : head;
                    if(newRow->nodeColors[other] == COLOR_FIRST){
                        numFirstColor++;
                    }else{
                        numSecondColor++;
                    }
                    iterEdge = getNextNodeEdge(dec,iterEdge,splitNode);
                }while(iterEdge != firstNodeEdge);

                COLOR_STATUS toNewNodeColor = numFirstColor < numSecondColor ? COLOR_FIRST : COLOR_SECOND;

                {

                    firstNodeEdge = getFirstNodeEdge(dec,splitNode);
                    iterEdge = firstNodeEdge;
                    do{
                        bool isCut = newRow->isEdgeCut[iterEdge];
                        spqr_node otherHead = findEdgeHead(dec, iterEdge);
                        spqr_node otherTail = findEdgeTail(dec, iterEdge);
                        spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
                        bool isMoveColor = newRow->nodeColors[otherEnd] == toNewNodeColor;
                        spqr_edge nextEdge = getNextNodeEdge(dec, iterEdge, splitNode); //Need to do this before we modify the edge :)

                        bool changeEdgeEnd = (isCut && isMoveColor) || (!isCut && !isMoveColor);
                        if(changeEdgeEnd){
                            if(otherHead == splitNode){
                                changeEdgeHead(dec,iterEdge,otherHead,newNode);
                            }else{
                                changeEdgeTail(dec,iterEdge,otherTail,newNode);
                            }
                        }
                        newRow->nodeColors[otherEnd] = UNCOLORED; //Clean up coloring information
                        //Ugly hack to make sure we can iterate neighbourhood whilst changing edge ends.
                        spqr_edge previousEdge = iterEdge;
                        iterEdge = nextEdge;
                        if(iterEdge == firstNodeEdge){
                            break;
                        }
                        if(changeEdgeEnd && previousEdge == firstNodeEdge){
                            firstNodeEdge = iterEdge;
                        }
                    }while(true);
                    newRow->reducedMembers[reducedMember].coloredNode = SPQR_INVALID_NODE;
                }

            }
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:{
            spqr_edge cutRepresentative = SPQR_INVALID_EDGE;
            MATREC_CALL(splitParallelMerging(dec,newRow,reducedMember,member,&newRowInfo->member,&cutRepresentative));
            spqr_node firstNode = SPQR_INVALID_NODE;
            spqr_node secondNode = SPQR_INVALID_NODE;
            spqr_node thirdNode = SPQR_INVALID_NODE;
            MATREC_CALL(createNode(dec,&firstNode));
            MATREC_CALL(createNode(dec,&secondNode));
            MATREC_CALL(createNode(dec,&thirdNode));
            spqr_edge first_edge = getFirstMemberEdge(dec, newRowInfo->member);
            spqr_edge edge = first_edge;

            do {
                if(edge != cutRepresentative){
                    setEdgeHeadAndTail(dec,edge,firstNode,secondNode);
                }else{
                    setEdgeHeadAndTail(dec,edge,firstNode,thirdNode);
                }
                edge = getNextMemberEdge(dec, edge);
            } while (edge != first_edge);

            newRowInfo->firstNode = secondNode;
            newRowInfo->secondNode = thirdNode;
            updateMemberType(dec, newRowInfo->member, SPQR_MEMBERTYPE_RIGID);

            break;
        }
        case SPQR_MEMBERTYPE_SERIES:{
            bool isCut = false;
            MATREC_CALL(splitSeriesMergingRowAddition(dec,newRow,reducedMember,member,&newRowInfo->member,&isCut));
            assert((newRow->reducedMembers[reducedMember].numChildren - newRow->reducedMembers[reducedMember].numPropagatedChildren) == 2);
            spqr_node firstNode = SPQR_INVALID_NODE;
            spqr_node secondNode = SPQR_INVALID_NODE;
            spqr_node thirdNode = SPQR_INVALID_NODE;
            spqr_node fourthNode = SPQR_INVALID_NODE;
            MATREC_CALL(createNode(dec,&firstNode));
            MATREC_CALL(createNode(dec,&secondNode));
            MATREC_CALL(createNode(dec,&thirdNode));
            MATREC_CALL(createNode(dec,&fourthNode));

            int reducedChildIndex = 0;

            spqr_edge reducedEdges[2];
            for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
                 i < newRow->reducedMembers[reducedMember].firstChild + newRow->reducedMembers[reducedMember].numChildren; ++i) {
                reduced_member_id child = newRow->childrenStorage[i];
                if(newRow->reducedMembers[child].type == TYPE_MERGED){
                    assert(reducedChildIndex < 2);
                    reducedEdges[reducedChildIndex] = markerOfParent(dec,newRow->reducedMembers[child].member);
                    reducedChildIndex++;
                }
            }

            spqr_edge first_edge = getFirstMemberEdge(dec, newRowInfo->member);
            spqr_edge edge = first_edge;


            do {

                if(edge == reducedEdges[0]){
                    setEdgeHeadAndTail(dec,edge,thirdNode,firstNode);
                }else if (edge == reducedEdges[1]){
                    if(isCut){
                        setEdgeHeadAndTail(dec,edge,fourthNode,secondNode);
                    }else{
                        setEdgeHeadAndTail(dec,edge,thirdNode,secondNode);
                    }
                }else{
                    setEdgeHeadAndTail(dec,edge,firstNode,secondNode);
                }
                edge = getNextMemberEdge(dec, edge);
            } while (edge != first_edge);

            newRowInfo->firstNode = thirdNode;
            newRowInfo->secondNode = fourthNode;

            updateMemberType(dec, newRowInfo->member, SPQR_MEMBERTYPE_RIGID);
            break;
        }
        default:
            assert(false);
            break;
    }
    return MATREC_OKAY;
}
static spqr_node getAdjacentNode(const NewRowInformation * const information,
                                 const spqr_node node1, const spqr_node node2){
    //Return the node which was passed which is one of the nodes in information
    if(node1 == information->firstNode || node1 == information->secondNode){
        return node1;
    }
    assert(node2 == information->firstNode || node2 == information->secondNode);
    return node2;
}

static MATREC_ERROR mergeRigidIntoParent(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                       const reduced_member_id reducedRigid, const spqr_member rigid,
                                       NewRowInformation * const information){

    //TODO: split up into smaller functions and re-use parts of series/parallel splitting
    spqr_member parent = information->member;
    spqr_edge parentToChild = markerOfParent(dec, rigid);
    spqr_edge childToParent = markerToParent(dec, rigid);

    spqr_node toParentHead = findEdgeHead(dec, childToParent);
    spqr_node toParentTail = findEdgeTail(dec, childToParent);
    assert(SPQRnodeIsValid(toParentHead) && SPQRnodeIsValid(toParentTail));
    assert(!NodePairHasTwo(&newRow->reducedMembers[reducedRigid].splitting_nodes) &&
           !NodePairIsEmpty(&newRow->reducedMembers[reducedRigid].splitting_nodes));
    spqr_node rigidSplit = newRow->reducedMembers[reducedRigid].splitting_nodes.first;
    spqr_node otherRigid = toParentHead == rigidSplit ? toParentTail : toParentHead;
    assert(otherRigid!= rigidSplit);

    spqr_node firstAdjacentNode = findEdgeHead(dec, parentToChild);
    spqr_node secondAdjacentNode = findEdgeTail(dec, parentToChild);
    assert(SPQRnodeIsValid(firstAdjacentNode) && SPQRnodeIsValid(secondAdjacentNode));
    spqr_node adjacentSplitNode = getAdjacentNode(information, firstAdjacentNode, secondAdjacentNode);
    spqr_node otherNode = adjacentSplitNode == firstAdjacentNode ? secondAdjacentNode : firstAdjacentNode;
    spqr_node otherSplitNode = adjacentSplitNode == information->firstNode ? information->secondNode : information->firstNode;
    assert(otherNode != information->firstNode && otherNode != information->secondNode);
    assert(otherSplitNode != firstAdjacentNode && otherSplitNode != secondAdjacentNode);

    spqr_node articulationEdgeHead = SPQR_INVALID_NODE;
#ifndef NDEBUG
    spqr_node articulationEdgeTail = SPQR_INVALID_NODE;
#endif
    if(newRow->reducedMembers[reducedRigid].allHaveCommonNode && (SPQRedgeIsValid(
            newRow->reducedMembers[reducedRigid].articulationEdge))){
        spqr_edge articulationEdge = newRow->reducedMembers[reducedRigid].articulationEdge;
        articulationEdgeHead = findEdgeHead(dec,articulationEdge);
#ifndef NDEBUG
        articulationEdgeTail = findEdgeTail(dec, articulationEdge);
#endif
        if(articulationEdge == childToParent ){
            swap_ints(&adjacentSplitNode,&otherSplitNode);
        }
    }

    COLOR_STATUS moveColor = newRow->nodeColors[otherRigid];
    //Remove the two marker edges which are merged (TODO: reuse edge memory?)
    {
        if(findEdgeHead(dec,childToParent) == rigidSplit){
            newRow->nodeColors[findEdgeTail(dec,childToParent)] = UNCOLORED;
        }else if (findEdgeTail(dec,childToParent) == rigidSplit){
            newRow->nodeColors[findEdgeHead(dec,childToParent)] = UNCOLORED;
        }
        removeEdgeFromMemberEdgeList(dec,childToParent,rigid);
        clearEdgeHeadAndTail(dec,childToParent);


        removeEdgeFromMemberEdgeList(dec,parentToChild,parent);
        clearEdgeHeadAndTail(dec,parentToChild); //TODO These functions call redundant finds

    }
    if(!newRow->reducedMembers[reducedRigid].allHaveCommonNode && newRow->reducedMembers[reducedRigid].numCutEdges > 0)  //articulation node splitting is easier to do before the merging
    {
        assert(moveColor == COLOR_FIRST || moveColor == COLOR_SECOND);
        //for each edge adjacent to the old rigid
        spqr_edge firstEdge = getFirstNodeEdge(dec, rigidSplit);
        spqr_edge iterEdge = firstEdge;
        do{
            bool isCut = newRow->isEdgeCut[iterEdge];
            spqr_node otherHead = findEdgeHead(dec, iterEdge);
            spqr_node otherTail = findEdgeTail(dec, iterEdge);
            spqr_node otherEnd = otherHead == rigidSplit ? otherTail : otherHead;
            bool isMoveColor = newRow->nodeColors[otherEnd] == moveColor;
            spqr_edge nextEdge = getNextNodeEdge(dec, iterEdge, rigidSplit); //Need to do this before we modify the edge :)

            bool changeEdgeEnd = (isCut && isMoveColor) || (!isCut && !isMoveColor);
            if(changeEdgeEnd){
                if(otherHead == rigidSplit){
                    changeEdgeHead(dec,iterEdge,otherHead,otherSplitNode);
                }else{
                    changeEdgeTail(dec,iterEdge,otherTail,otherSplitNode);
                }
            }
            newRow->nodeColors[otherEnd] = UNCOLORED;

            //Ugly hack to make sure we can iterate neighbourhood whilst changing edge ends.
            spqr_edge previousEdge = iterEdge;
            iterEdge = nextEdge;
            if(iterEdge == firstEdge){
                break;
            }
            if(changeEdgeEnd && previousEdge == firstEdge){
                firstEdge = iterEdge;
            }
        }while(true);
        newRow->reducedMembers[reducedRigid].coloredNode = SPQR_INVALID_NODE;
    }

    //Identify the members with each other
    {
        spqr_member newMember = mergeMembers(dec, rigid, parent);
        spqr_member toRemoveFrom = newMember == rigid ? parent : rigid;

        mergeMemberEdgeList(dec, newMember, toRemoveFrom);
        if (toRemoveFrom == parent) { //Correctly update the parent information
            updateMemberParentInformation(dec,newMember,toRemoveFrom);
        }
        updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);
        information->member = newMember;
    }

    //identify rigid_split with adjacent_split and other_rigid with other_node
    spqr_node mergedSplit = mergeNodes(dec, rigidSplit, adjacentSplitNode);
    spqr_node splitToRemove = mergedSplit == rigidSplit ? adjacentSplitNode : rigidSplit;
    if(splitToRemove == information->firstNode){
        information->firstNode = mergedSplit;
    }else if (splitToRemove == information->secondNode){
        information->secondNode = mergedSplit;
    }
    mergeNodes(dec,otherRigid,otherNode); //Returns the update node ID, but we do not need this here.

    if(newRow->reducedMembers[reducedRigid].allHaveCommonNode){
        if(SPQRedgeIsValid(newRow->reducedMembers[reducedRigid].articulationEdge)) {
            if(newRow->reducedMembers[reducedRigid].articulationEdge != childToParent) {
                spqr_edge articulationEdge = newRow->reducedMembers[reducedRigid].articulationEdge;
                if (articulationEdgeHead == rigidSplit) {
                    changeEdgeHead(dec, articulationEdge, mergedSplit, otherSplitNode);
                } else {
                    assert(articulationEdgeTail == rigidSplit);
                    changeEdgeTail(dec, articulationEdge, mergedSplit, otherSplitNode);
                }
            }
        }else{
            cut_edge_id cutEdgeIdx = newRow->reducedMembers[reducedRigid].firstCutEdge;
            do {
                spqr_edge cutEdge = newRow->cutEdges[cutEdgeIdx].edge;
                spqr_node edgeHead = findEdgeHead(dec, cutEdge);
                if (edgeHead == mergedSplit) {
                    changeEdgeHead(dec, cutEdge, edgeHead, otherSplitNode);
                } else {
                    changeEdgeTail(dec, cutEdge, findEdgeTail(dec, cutEdge), otherSplitNode);
                }

                cutEdgeIdx = newRow->cutEdges[cutEdgeIdx].nextMember;
            } while (cutEdgeIsValid(cutEdgeIdx));
        }
        return MATREC_OKAY;
    }
    return MATREC_OKAY;
}

static spqr_member mergeParallelIntoParent(MATRECGraphicDecomposition * dec,
                                           const spqr_member member, const spqr_member parent,
                                           NewRowInformation * const information,
                                           const spqr_edge cutRepresentative){

    spqr_edge parentToChild = markerOfParent(dec, member);
    assert(findEdgeMemberNoCompression(dec,parentToChild) == parent);

    //Remove the marker edges which are merged
    //TODO: reuse edge memory?
    removeEdgeFromMemberEdgeList(dec,parentToChild,parent);

    spqr_edge childToParent = markerToParent(dec, member);
    removeEdgeFromMemberEdgeList(dec,childToParent,member);

    {
        spqr_node firstAdjacentNode = findEdgeHead(dec, parentToChild);
        spqr_node secondAdjacentNode = findEdgeTail(dec, parentToChild);
        assert(SPQRnodeIsValid(firstAdjacentNode) && SPQRnodeIsValid(secondAdjacentNode));
        clearEdgeHeadAndTail(dec,parentToChild); //By the merging procedure,the parent is always a rigid member.

        spqr_node adjacentSplitNode = getAdjacentNode(information, firstAdjacentNode, secondAdjacentNode);
        spqr_node otherNode = adjacentSplitNode == firstAdjacentNode ? secondAdjacentNode : firstAdjacentNode;
        spqr_node otherSplitNode = adjacentSplitNode == information->firstNode ? information->secondNode : information->firstNode;

        assert(otherNode != information->firstNode && otherNode != information->secondNode);
        assert(otherSplitNode != firstAdjacentNode && otherSplitNode != secondAdjacentNode);

        spqr_edge first_edge = getFirstMemberEdge(dec, member);
        spqr_edge edge = first_edge;

        do {
            if(edge == cutRepresentative){
                setEdgeHeadAndTail(dec,edge,otherSplitNode,otherNode);
            }else{
                setEdgeHeadAndTail(dec,edge,adjacentSplitNode,otherNode);
            }
            edge = getNextMemberEdge(dec, edge);
        } while (edge != first_edge);

    }

    spqr_member newMember = mergeMembers(dec, member, parent);
    spqr_member toRemoveFrom = newMember == member ? parent : member;

    mergeMemberEdgeList(dec,newMember,toRemoveFrom);
    if(toRemoveFrom == parent){ //Correctly update the parent information
        updateMemberParentInformation(dec,newMember,toRemoveFrom);
    }
    updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);

    return newMember;
}

static spqr_edge seriesGetOtherEdge(const MATRECGraphicDecomposition * dec, const MATRECGraphicRowAddition * newRow,
                                    const reduced_member_id reducedMember){
    for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
         i < newRow->reducedMembers[reducedMember].firstChild +
             newRow->reducedMembers[reducedMember].numChildren; ++i) {
        reduced_member_id child = newRow->childrenStorage[i];
        if (newRow->reducedMembers[child].type == TYPE_MERGED) {
            return markerOfParent(dec, newRow->reducedMembers[child].member);
        }
    }
    assert(false);
    return SPQR_INVALID_EDGE;
}
static MATREC_ERROR mergeSeriesIntoParent(MATRECGraphicDecomposition * dec,
                                        const spqr_member member, spqr_member * const parent,
                                        NewRowInformation * const information, const bool isCut, const spqr_edge childDecompositionEdge){
    assert(getNumMemberEdges(dec,member) == 3);
    spqr_edge parentToChild = markerOfParent(dec, member);
    assert(findEdgeMemberNoCompression(dec,parentToChild) == *parent);

    //Remove the marker edges which are merged
    removeEdgeFromMemberEdgeList(dec,parentToChild,*parent);
    spqr_edge childToParent = markerToParent(dec, member);
    removeEdgeFromMemberEdgeList(dec,childToParent,member);

    {
        spqr_node newNode = SPQR_INVALID_NODE;
        MATREC_CALL(createNode(dec,&newNode));

        spqr_node firstAdjacentNode = findEdgeHead(dec, parentToChild);
        spqr_node secondAdjacentNode = findEdgeTail(dec, parentToChild);
        assert(SPQRnodeIsValid(firstAdjacentNode) && SPQRnodeIsValid(secondAdjacentNode));
        clearEdgeHeadAndTail(dec,parentToChild);

        spqr_node adjacentSplitNode = getAdjacentNode(information, firstAdjacentNode, secondAdjacentNode);
        spqr_node otherNode = adjacentSplitNode == firstAdjacentNode ? secondAdjacentNode : firstAdjacentNode;
        spqr_node otherSplitNode = adjacentSplitNode == information->firstNode ? information->secondNode : information->firstNode;

        assert(otherNode != information->firstNode && otherNode != information->secondNode);
        assert(otherSplitNode != firstAdjacentNode && otherSplitNode != secondAdjacentNode);

        spqr_edge first_edge = getFirstMemberEdge(dec, member);
        spqr_edge edge = first_edge;

        do {
            if(edge == childDecompositionEdge){
                if(isCut){
                    setEdgeHeadAndTail(dec,edge,otherSplitNode,newNode);
                }else{
                    setEdgeHeadAndTail(dec,edge,adjacentSplitNode,newNode);
                }
            }else{
                setEdgeHeadAndTail(dec,edge,otherNode,newNode);
            }

            edge = getNextMemberEdge(dec, edge);
        } while (edge != first_edge);

    }

    spqr_member newMember = mergeMembers(dec, member, *parent);
    spqr_member toRemoveFrom = newMember == member ? *parent : member;

    mergeMemberEdgeList(dec,newMember,toRemoveFrom);
    if(toRemoveFrom == *parent){ //Correctly update the parent information
        updateMemberParentInformation(dec,newMember,toRemoveFrom);
    }
    updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);
    *parent = newMember;
    return MATREC_OKAY;
}

static MATREC_ERROR mergeIntoLargeComponent(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition * newRow,
                                          const reduced_member_id reducedMember, NewRowInformation * const newRowInformation){
    spqr_member member = newRow->reducedMembers[reducedMember].member;
    assert(findMemberParentNoCompression(dec,member) == newRowInformation->member);
    spqr_member memberForMerging = member;

    switch(getMemberType(dec,member)){
        case SPQR_MEMBERTYPE_RIGID:{
            MATREC_CALL(mergeRigidIntoParent(dec, newRow, reducedMember, member, newRowInformation));
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:{
            spqr_edge cutRepresentative = SPQR_INVALID_EDGE;
            MATREC_CALL(splitParallelMerging(dec,newRow,reducedMember,member,&memberForMerging,&cutRepresentative));
            assert(findMemberParentNoCompression(dec,memberForMerging) == newRowInformation->member);
            spqr_member newID = mergeParallelIntoParent(dec, memberForMerging, newRowInformation->member, newRowInformation, cutRepresentative);
            newRowInformation->member = newID;
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:
        {
            bool isCut = false;
            MATREC_CALL(splitSeriesMergingRowAddition(dec,newRow,reducedMember,member,&memberForMerging,&isCut));
            assert(findMemberParentNoCompression(dec,memberForMerging) == newRowInformation->member);

            spqr_edge otherRepresentative = seriesGetOtherEdge(dec, newRow, reducedMember);
            MATREC_CALL(mergeSeriesIntoParent(dec,memberForMerging,&newRowInformation->member,newRowInformation, isCut,otherRepresentative));
            break;
        }

        default:
            break;
    }
    return MATREC_OKAY;
}

static MATREC_ERROR mergeTree(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow, reduced_member_id root,
                            NewRowInformation * const newRowInformation){

    //Multiple members: we need to merge them together
    //We start by transforming the root, and then iteratively merging any of the relevant children into it
    MATREC_CALL(transformRoot(dec, newRow, root,newRowInformation));

    //Iteratively merge into root component using DFS

    //We reuse the data for determining the types, which has similar call stack data and uses more memories
    assert(newRow->memMergeTreeCallData >= newRow->numReducedMembers);

    int depth = 0;
    MergeTreeCallData * stack = newRow->mergeTreeCallData;

    stack[0].currentChild = newRow->reducedMembers[root].firstChild;
    stack[0].id = root;
    do{
        if(stack[depth].currentChild == newRow->reducedMembers[stack[depth].id].firstChild +
                                        newRow->reducedMembers[stack[depth].id].numChildren){
            --depth;
            continue;
        }
        reduced_member_id reducedChild = newRow->childrenStorage[stack[depth].currentChild];
        ++stack[depth].currentChild;
        if(newRow->reducedMembers[reducedChild].type == TYPE_MERGED){
            MATREC_CALL(mergeIntoLargeComponent(dec,newRow, reducedChild, newRowInformation));
            ++depth;
            assert(depth < newRow->memMergeTreeCallData);
            stack[depth].id = reducedChild;
            stack[depth].currentChild = newRow->reducedMembers[reducedChild].firstChild;
        }

    }while(depth >= 0);
    return MATREC_OKAY;
}
static MATREC_ERROR transformComponentRowAddition(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow,
                                     MATRECRowReducedComponent* component, NewRowInformation * const newRowInformation){
    assert(component);
    if(newRow->reducedMembers[component->root].numChildren == newRow->reducedMembers[component->root].numPropagatedChildren){
        //No merging necessary, only a single component
        reduced_member_id reducedMember = component->root;
        assert(reducedMemberIsValid(reducedMember));
        spqr_member member = newRow->reducedMembers[reducedMember].member;
        SPQRMemberType type = getMemberType(dec, member);

        switch(type){
            case SPQR_MEMBERTYPE_RIGID:
                MATREC_CALL(transformSingleRigid(dec,newRow,reducedMember,member,newRowInformation));
                break;
            case SPQR_MEMBERTYPE_PARALLEL: {
                MATREC_CALL(transformSingleParallel(dec,newRow,reducedMember,member,&newRowInformation->member));
                break;
            }
            case SPQR_MEMBERTYPE_SERIES: {
                newRowInformation->member = member;
                break;
            }
            case SPQR_MEMBERTYPE_LOOP:{
                int numEdges = getNumMemberEdges(dec,member);
                if( numEdges == 2){
                    changeLoopToSeries(dec,member);
                }
                newRowInformation->member = member;
                break;
            }
            default:
                assert(false);
                break;
        }

        return MATREC_OKAY;
    }

    MATREC_CALL(mergeTree(dec,newRow,component->root,newRowInformation));

    return MATREC_OKAY;
}


MATREC_ERROR MATRECcreateGraphicRowAddition(MATREC* env, MATRECGraphicRowAddition** pNewRow ){
    assert(env);
    MATREC_CALL(MATRECallocBlock(env,pNewRow));
    MATRECGraphicRowAddition * newRow = *pNewRow;

    newRow->remainsGraphic = true;

    newRow->reducedMembers = NULL;
    newRow->memReducedMembers = 0;
    newRow->numReducedMembers = 0;

    newRow->reducedComponents = NULL;
    newRow->memReducedComponents = 0;
    newRow->numReducedComponents = 0;

    newRow->memberInformation = NULL;
    newRow->memMemberInformation = 0;
    newRow->numMemberInformation = 0;

    newRow->cutEdges = NULL;
    newRow->memCutEdges = 0;
    newRow->numCutEdges = 0;
    newRow->firstOverallCutEdge = INVALID_CUT_EDGE;

    newRow->childrenStorage = NULL;
    newRow->memChildrenStorage = 0;
    newRow->numChildrenStorage = 0;

    newRow->newRowIndex = MATREC_INVALID_ROW;

    newRow->newColumnEdges = NULL;
    newRow->memColumnEdges = 0;
    newRow->numColumnEdges = 0;

    newRow->leafMembers = NULL;
    newRow->numLeafMembers = 0;
    newRow->memLeafMembers = 0;

    newRow->decompositionColumnEdges = NULL;
    newRow->memDecompositionColumnEdges = 0;
    newRow->numDecompositionColumnEdges = 0;

    newRow->isEdgeCut = NULL;
    newRow->memIsEdgeCut = 0;
    newRow->numIsEdgeCut = 0;

    newRow->nodeColors = NULL;
    newRow->memNodeColors = 0;

    newRow->articulationNodes = NULL;
    newRow->memArticulationNodes = 0;
    newRow->numArticulationNodes = 0;

    newRow->articulationNodeSearchInfo = NULL;
    newRow->memNodeSearchInfo = 0;

    newRow->crossingPathCount = NULL;
    newRow->memCrossingPathCount = 0;

    newRow->intersectionDFSData = NULL;
    newRow->memIntersectionDFSData = 0;

    newRow->colorDFSData = NULL;
    newRow->memColorDFSData = 0;

    newRow->artDFSData = NULL;
    newRow->memArtDFSData = 0;

    newRow->createReducedMembersCallstack = NULL;
    newRow->memCreateReducedMembersCallstack = 0;

    newRow->intersectionPathDepth = NULL;
    newRow->memIntersectionPathDepth = 0;

    newRow->intersectionPathParent = NULL;
    newRow->memIntersectionPathParent = 0;

    newRow->mergeTreeCallData = NULL;
    newRow->memMergeTreeCallData = 0;

    return MATREC_OKAY;
}

void MATRECfreeGraphicRowAddition(MATREC* env, MATRECGraphicRowAddition ** pNewRow){
    assert(*pNewRow);

    MATRECGraphicRowAddition * newRow = *pNewRow;
    //TODO: check if everything is truly freed in reverse order

    MATRECfreeBlockArray(env,&newRow->createReducedMembersCallstack);
    MATRECfreeBlockArray(env,&newRow->artDFSData);
    MATRECfreeBlockArray(env,&newRow->colorDFSData);
    MATRECfreeBlockArray(env,&newRow->mergeTreeCallData);
    MATRECfreeBlockArray(env,&newRow->intersectionDFSData);
    MATRECfreeBlockArray(env,&newRow->intersectionPathParent);
    MATRECfreeBlockArray(env,&newRow->intersectionPathDepth);
    MATRECfreeBlockArray(env,&newRow->crossingPathCount);
    MATRECfreeBlockArray(env,&newRow->articulationNodeSearchInfo);
    MATRECfreeBlockArray(env,&newRow->articulationNodes);
    MATRECfreeBlockArray(env,&newRow->nodeColors);
    MATRECfreeBlockArray(env,&newRow->isEdgeCut);
    if(newRow->decompositionColumnEdges){
        MATRECfreeBlockArray(env,&newRow->decompositionColumnEdges);
    }
    if(newRow->newColumnEdges){
        MATRECfreeBlockArray(env,&newRow->newColumnEdges);
    }
    if(newRow->childrenStorage){
        MATRECfreeBlockArray(env,&newRow->childrenStorage);
    }
    if(newRow->cutEdges){
        MATRECfreeBlockArray(env,&newRow->cutEdges);
    }
    if(newRow->memberInformation){
        MATRECfreeBlockArray(env,&newRow->memberInformation);
    }
    if(newRow->reducedComponents){
        MATRECfreeBlockArray(env,&newRow->reducedComponents);
    }
    if(newRow->reducedMembers){
        MATRECfreeBlockArray(env,&newRow->reducedMembers);
    }
    MATRECfreeBlockArray(env,&newRow->leafMembers);
    MATRECfreeBlock(env,pNewRow);
}

MATREC_ERROR MATRECGraphicRowAdditionCheck(MATRECGraphicDecomposition * dec, MATRECGraphicRowAddition * newRow, const MATREC_row row, const MATREC_col * columns, size_t numColumns){
    assert(dec);
    assert(newRow);
    assert(numColumns == 0 || columns );

    newRow->remainsGraphic = true;
    cleanUpPreviousIteration(dec,newRow);

    MATREC_CALL(newRowUpdateRowInformation(dec,newRow,row,columns,numColumns));
    MATREC_CALL(constructRowReducedDecomposition(dec,newRow));
    MATREC_CALL(createReducedDecompositionCutEdges(dec,newRow));

    MATREC_CALL(determineLeafReducedMembers(dec,newRow));
    MATREC_CALL(allocateRigidSearchMemory(dec,newRow));
    MATREC_CALL(allocateTreeSearchMemory(dec,newRow));
    //Check for each component if the cut edges propagate through a row tree marker to a cut edge in another component
    //From the leafs inward.
    propagateComponents(dec,newRow);
    //It can happen that we are not graphic by some of the checked components.
    //In that case, further checking may lead to errors as some invariants that the code assumes will be broken.
    if(newRow->remainsGraphic){
        for (int i = 0; i < newRow->numReducedComponents; ++i) {
            determineMergeableTypes(dec,newRow,newRow->reducedComponents[i].root);
            //exit early if one is not graphic
            if(!newRow->remainsGraphic){
                break;
            }
        }
    }

    cleanUpRowMemberInformation(newRow);
    if(!newRow->remainsGraphic){
        for (int i = 0; i < newRow->numReducedMembers; ++i) {
            if(SPQRnodeIsValid(newRow->reducedMembers[i].coloredNode)){
                zeroOutColors(dec,newRow,newRow->reducedMembers[i].coloredNode);
                newRow->reducedMembers[i].coloredNode = SPQR_INVALID_NODE;

            }
        }
    }
    return MATREC_OKAY;
}

MATREC_ERROR MATRECGraphicRowAdditionAdd(MATRECGraphicDecomposition *dec, MATRECGraphicRowAddition *newRow){
    assert(newRow->remainsGraphic);
    if(newRow->numReducedComponents == 0){
        spqr_member newMember = SPQR_INVALID_MEMBER;
        MATREC_CALL(createStandaloneParallel(dec,newRow->newColumnEdges,newRow->numColumnEdges,newRow->newRowIndex,&newMember));
    }else if (newRow->numReducedComponents == 1){
        NewRowInformation information = emptyNewRowInformation();
        MATREC_CALL(transformComponentRowAddition(dec,newRow,&newRow->reducedComponents[0],&information));

        if(newRow->numColumnEdges == 0){
            spqr_edge row_edge = SPQR_INVALID_EDGE;
            MATREC_CALL(createRowEdge(dec,information.member,&row_edge,newRow->newRowIndex));
            if(SPQRnodeIsValid(information.firstNode)){
                setEdgeHeadAndTail(dec,row_edge,information.firstNode,information.secondNode);
            }
        }else{
            spqr_member new_row_parallel = SPQR_INVALID_MEMBER;
            MATREC_CALL(createConnectedParallel(dec,newRow->newColumnEdges,newRow->numColumnEdges,newRow->newRowIndex,&new_row_parallel));
            spqr_edge markerEdge = SPQR_INVALID_EDGE;
            spqr_edge ignore = SPQR_INVALID_EDGE;
            MATREC_CALL(createMarkerPairWithReferences(dec,information.member,new_row_parallel,true,&markerEdge,&ignore));
            if(SPQRnodeIsValid(information.firstNode)){
                setEdgeHeadAndTail(dec,markerEdge,information.firstNode,information.secondNode);
            }
        }
        if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
            assert(getNumMemberEdges(dec,information.member) == 2 || getNumMemberEdges(dec,information.member) == 3);
            if(getNumMemberEdges(dec,information.member) == 3){
                changeLoopToSeries(dec,information.member);
            }
        }
    }else{

#ifndef NDEBUG
        int numDecComponentsBefore = numConnectedComponents(dec);
#endif
        spqr_member new_row_parallel = SPQR_INVALID_MEMBER;
        MATREC_CALL(createConnectedParallel(dec,newRow->newColumnEdges,newRow->numColumnEdges,newRow->newRowIndex,&new_row_parallel));
        for (int i = 0; i < newRow->numReducedComponents; ++i) {
            NewRowInformation information = emptyNewRowInformation();
            MATREC_CALL(transformComponentRowAddition(dec,newRow,&newRow->reducedComponents[i],&information));
            if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
                assert(getNumMemberEdges(dec,information.member) == 1);
                moveEdgeToNewMember(dec, getFirstMemberEdge(dec,information.member),information.member,new_row_parallel);
                dec->members[information.member].type = SPQR_MEMBERTYPE_UNASSIGNED;
            }else{
                reorderComponent(dec,information.member); //Make sure the new component is the root of the local decomposition tree
                spqr_edge markerEdge = SPQR_INVALID_EDGE;
                spqr_edge ignore = SPQR_INVALID_EDGE;
                MATREC_CALL(createMarkerPairWithReferences(dec,new_row_parallel,information.member,false,&ignore,&markerEdge));
                if(SPQRnodeIsValid(information.firstNode)){
                    setEdgeHeadAndTail(dec,markerEdge,information.firstNode,information.secondNode);
                }
            }

        }
        decreaseNumConnectedComponents(dec,newRow->numReducedComponents-1);
        assert(numConnectedComponents(dec) == (numDecComponentsBefore - newRow->numReducedComponents + 1));
    }
    for (int i = 0; i < newRow->numReducedMembers; ++i) {
        if(SPQRnodeIsValid(newRow->reducedMembers[i].coloredNode)){
            zeroOutColors(dec,newRow,newRow->reducedMembers[i].coloredNode);
            newRow->reducedMembers[i].coloredNode = SPQR_INVALID_NODE;

        }
    }
    return MATREC_OKAY;
}

bool MATRECGraphicRowAdditionRemainsGraphic(const MATRECGraphicRowAddition *newRow){
    return newRow->remainsGraphic;
}
