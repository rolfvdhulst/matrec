#include "matrec/Network.h"
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

typedef int spqr_arc;
#define SPQR_INVALID_ARC (-1)

static bool SPQRarcIsInvalid(spqr_arc arc){
    return arc <  0;
}
static bool SPQRarcIsValid(spqr_arc arc){
    return !SPQRarcIsInvalid(arc);
}

typedef enum {
    SPQR_MEMBERTYPE_RIGID = 0, //Also known as triconnected components
    SPQR_MEMBERTYPE_PARALLEL = 1,//Also known as a 'bond'
    SPQR_MEMBERTYPE_SERIES = 2, //Also known as 'polygon' or 'cycle'
    SPQR_MEMBERTYPE_LOOP = 3,
    SPQR_MEMBERTYPE_UNASSIGNED = 4 // To indicate that the member has been merged/is not representative; this is just there to catch errors.
} SPQRMemberType;

typedef struct{
    spqr_arc previous;
    spqr_arc next;
}MATRECNetworkDecompositionArcListNode;

typedef struct {
    spqr_node representativeNode;
    spqr_arc firstArc;//first arc of the neighbouring arcs
    int numArcs;
} MATRECNetworkDecompositionNode;

typedef struct {
    spqr_node head;
    spqr_node tail;
    spqr_member member;
    spqr_member childMember;
    MATRECNetworkDecompositionArcListNode headArcListNode;
    MATRECNetworkDecompositionArcListNode tailArcListNode;
    MATRECNetworkDecompositionArcListNode arcListNode; //Linked-list node of the array of arcs of the member which this arc is in

    spqr_element element;

    //Signed union-find for arc directions
    //For non-rigid members every arc is it's own representative, and the direction is simply given by the boolean
    //For rigid members, every arc is represented by another arc in the member,
    //and the direction can be found by multiplying the signs along the union-find path
    spqr_arc representative;
    bool reversed;
} MATRECNetworkDecompositionArc;

typedef struct {
    spqr_member representativeMember;
    SPQRMemberType type;

    spqr_member parentMember;
    spqr_arc markerToParent;
    spqr_arc markerOfParent;

    spqr_arc firstArc; //First of the members' linked-list arc array
    int numArcs;
} MATRECNetworkDecompositionMember;

struct MATRECNetworkDecompositionImpl {
    int numArcs;
    int memArcs;
    MATRECNetworkDecompositionArc *arcs;
    spqr_arc firstFreeArc;

    int memMembers;
    int numMembers;
    MATRECNetworkDecompositionMember *members;

    int memNodes;
    int numNodes;
    MATRECNetworkDecompositionNode *nodes;

    int memRows;
    int numRows;
    spqr_arc * rowArcs;

    int memColumns;
    int numColumns;
    spqr_arc * columnArcs;

    MATREC * env;

    int numConnectedComponents;
};

static void swap_ints(int* a, int* b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

static bool nodeIsRepresentative(const MATRECNetworkDecomposition *dec, spqr_node node) {
    assert(dec);
    assert(node < dec->memNodes);
    assert(SPQRnodeIsValid(node));

    return SPQRnodeIsInvalid(dec->nodes[node].representativeNode);
}

static spqr_node findNode(MATRECNetworkDecomposition *dec, spqr_node node) {
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

static spqr_node findNodeNoCompression(const MATRECNetworkDecomposition *dec, spqr_node node) {
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

static spqr_node findArcTail(MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_node representative = findNode(dec, dec->arcs[arc].tail);
    dec->arcs[arc].tail = representative; //update the arc information

    return representative;
}
static spqr_node findArcHead(MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_node representative = findNode(dec, dec->arcs[arc].head);
    dec->arcs[arc].head = representative;//update the arc information

    return representative;
}
static spqr_node findArcHeadNoCompression(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_node representative = findNodeNoCompression(dec, dec->arcs[arc].head);
    return representative;
}
static spqr_node findArcTailNoCompression(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_node representative = findNodeNoCompression(dec, dec->arcs[arc].tail);
    return representative;
}


static spqr_arc getFirstNodeArc(const MATRECNetworkDecomposition * dec, spqr_node node){
    assert(dec);
    assert(SPQRnodeIsValid(node));
    assert(node < dec->memNodes);
    return dec->nodes[node].firstArc;
}
static spqr_arc getNextNodeArcNoCompression(const MATRECNetworkDecomposition * dec, spqr_arc arc, spqr_node node){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    assert(nodeIsRepresentative(dec,node));

    if(findArcHeadNoCompression(dec,arc) == node){
        arc = dec->arcs[arc].headArcListNode.next;
    }else{
        assert(findArcTailNoCompression(dec,arc) == node);
        arc = dec->arcs[arc].tailArcListNode.next;
    }
    return arc;
}
static spqr_arc getNextNodeArc(MATRECNetworkDecomposition * dec, spqr_arc arc, spqr_node node){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    assert(nodeIsRepresentative(dec,node));

    if(findArcHead(dec,arc) == node){
        arc = dec->arcs[arc].headArcListNode.next;
    }else{
        assert(findArcTailNoCompression(dec,arc) == node);
        dec->arcs[arc].tail = node; //This assignment is not necessary but speeds up future queries.
        arc = dec->arcs[arc].tailArcListNode.next;
    }
    return arc;
}
static spqr_arc getPreviousNodeArc(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_node node){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    assert(nodeIsRepresentative(dec,node));

    if(findArcHead(dec,arc) == node){
        arc = dec->arcs[arc].headArcListNode.previous;
    }else{
        assert(findArcTailNoCompression(dec,arc) == node);
        dec->arcs[arc].tail = node; //This assignment is not necessary but speeds up future queries.
        arc = dec->arcs[arc].tailArcListNode.previous;
    }
    return arc;
}

static void mergeNodeArcList(MATRECNetworkDecomposition *dec, spqr_node toMergeInto, spqr_node toRemove){

    spqr_arc firstIntoArc = getFirstNodeArc(dec, toMergeInto);
    spqr_arc firstFromArc = getFirstNodeArc(dec, toRemove);
    if(SPQRarcIsInvalid(firstIntoArc)){
        //new node has no arcs
        dec->nodes[toMergeInto].numArcs += dec->nodes[toRemove].numArcs;
        dec->nodes[toRemove].numArcs = 0;

        dec->nodes[toMergeInto].firstArc = dec->nodes[toRemove].firstArc;
        dec->nodes[toRemove].firstArc = SPQR_INVALID_ARC;

        return;
    }else if (SPQRarcIsInvalid(firstFromArc)){
        //Old node has no arcs; we can just return
        return;
    }

    spqr_arc lastIntoArc = getPreviousNodeArc(dec, firstIntoArc, toMergeInto);
    assert(SPQRarcIsValid(lastIntoArc));
    spqr_arc lastFromArc = getPreviousNodeArc(dec, firstFromArc, toRemove);
    assert(SPQRarcIsValid(lastFromArc));


    MATRECNetworkDecompositionArcListNode * firstIntoNode = findArcHead(dec, firstIntoArc) == toMergeInto ?
                                                           &dec->arcs[firstIntoArc].headArcListNode :
                                                           &dec->arcs[firstIntoArc].tailArcListNode;
    MATRECNetworkDecompositionArcListNode * lastIntoNode = findArcHead(dec, lastIntoArc) == toMergeInto ?
                                                          &dec->arcs[lastIntoArc].headArcListNode :
                                                          &dec->arcs[lastIntoArc].tailArcListNode;

    MATRECNetworkDecompositionArcListNode * firstFromNode = findArcHead(dec, firstFromArc) == toRemove ?
                                                           &dec->arcs[firstFromArc].headArcListNode :
                                                           &dec->arcs[firstFromArc].tailArcListNode;
    MATRECNetworkDecompositionArcListNode * lastFromNode = findArcHead(dec, lastFromArc) == toRemove ?
                                                          &dec->arcs[lastFromArc].headArcListNode :
                                                          &dec->arcs[lastFromArc].tailArcListNode;

    firstIntoNode->previous = lastFromArc;
    lastIntoNode->next = firstFromArc;
    firstFromNode->previous = lastIntoArc;
    lastFromNode->next = firstIntoArc;

    dec->nodes[toMergeInto].numArcs += dec->nodes[toRemove].numArcs;
    dec->nodes[toRemove].numArcs = 0;
    dec->nodes[toRemove].firstArc = SPQR_INVALID_ARC;
}

static void arcFlipReversed(MATRECNetworkDecomposition *dec, spqr_arc arc){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    dec->arcs[arc].reversed = !dec->arcs[arc].reversed;
}
static void arcSetReversed(MATRECNetworkDecomposition *dec, spqr_arc arc, bool reversed){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    dec->arcs[arc].reversed = reversed;
}
static void arcSetRepresentative(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_arc representative){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    assert(representative == SPQR_INVALID_ARC || SPQRarcIsValid(representative));
    dec->arcs[arc].representative = representative;
}

static spqr_node mergeNodes(MATRECNetworkDecomposition *dec, spqr_node first, spqr_node second) {
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
    //first becomes representative; we merge all of the arcs of second into first
    mergeNodeArcList(dec,first,second);
    dec->nodes[second].representativeNode = first;
    if (firstRank == secondRank) {
        --dec->nodes[first].representativeNode;
    }
    return first;
}

static bool memberIsRepresentative(const MATRECNetworkDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(member < dec->memMembers);
    assert(SPQRmemberIsValid(member));

    return SPQRmemberIsInvalid(dec->members[member].representativeMember);
}

static spqr_member findMember(MATRECNetworkDecomposition *dec, spqr_member member) {
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

static spqr_member findMemberNoCompression(const MATRECNetworkDecomposition *dec, spqr_member member) {
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

static spqr_member mergeMembers(MATRECNetworkDecomposition *dec, spqr_member first, spqr_member second) {
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

static spqr_member findArcMember(MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_member representative = findMember(dec, dec->arcs[arc].member);
    dec->arcs[arc].member = representative;
    return representative;
}

static spqr_member findArcMemberNoCompression(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_member representative = findMemberNoCompression(dec, dec->arcs[arc].member);
    return representative;
}

static spqr_member findMemberParent(MATRECNetworkDecomposition *dec, spqr_member member) {
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

static spqr_member findMemberParentNoCompression(const MATRECNetworkDecomposition *dec, spqr_member member) {
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

static spqr_member findArcChildMember(MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_member representative = findMember(dec, dec->arcs[arc].childMember);
    dec->arcs[arc].childMember = representative;
    return representative;
}

static spqr_member findArcChildMemberNoCompression(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    spqr_member representative = findMemberNoCompression(dec, dec->arcs[arc].childMember);
    return representative;
}

//TODO: fix usages, is misleading. Only accounts for CHILD markers, not parent markers!
static bool arcIsMarker(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    return SPQRmemberIsValid(dec->arcs[arc].childMember);
}

static bool arcIsTree(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    return SPQRelementIsRow(dec->arcs[arc].element);
}

typedef struct {
    spqr_arc representative;
    bool reversed;
} ArcSign;
//find
static bool arcIsRepresentative(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(arc < dec->memArcs);
    assert(SPQRarcIsValid(arc));

    return SPQRarcIsInvalid(dec->arcs[arc].representative);
}

static ArcSign findArcSign(MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(arc < dec->memArcs);
    assert(SPQRarcIsValid(arc));

    spqr_arc current = arc;
    spqr_arc next;

    bool totalReversed = dec->arcs[current].reversed;
    //traverse down tree to find the root
    while (SPQRarcIsValid(next = dec->arcs[current].representative)) {
        current = next;
        assert(current < dec->memArcs);
        //swap boolean only if new arc is reversed
        totalReversed = (totalReversed != dec->arcs[current].reversed);
    }

    spqr_arc root = current;
    current = arc;

    bool currentReversed = totalReversed != dec->arcs[root].reversed;
    //update all pointers along path to point to root, flattening the tree

    while (SPQRarcIsValid(next = dec->arcs[current].representative)) {
        bool wasReversed = dec->arcs[current].reversed;

        dec->arcs[current].reversed = currentReversed;
        currentReversed = (currentReversed != wasReversed);

        dec->arcs[current].representative = root;
        current = next;
        assert(current < dec->memArcs);
    }

    ArcSign sign;
    sign.reversed = totalReversed;
    sign.representative = root;
    return sign;
}

static ArcSign findArcSignNoCompression(const MATRECNetworkDecomposition *dec, spqr_arc arc) {
    assert(dec);
    assert(arc < dec->memArcs);
    assert(SPQRarcIsValid(arc));

    spqr_arc current = arc;
    spqr_arc next;

    bool totalReversed = dec->arcs[current].reversed;
    //traverse down tree to find the root
    while (SPQRarcIsValid(next = dec->arcs[current].representative)) {
        current = next;
        assert(current < dec->memArcs);
        //swap boolean only if new arc is reversed
        totalReversed = (totalReversed != dec->arcs[current].reversed);
    }
    ArcSign sign;
    sign.reversed = totalReversed;
    sign.representative = current;
    return sign;
}
//Find the arc tail/head, but accounting for reflection
static spqr_node findEffectiveArcHead(MATRECNetworkDecomposition *dec, spqr_arc arc){
    assert(dec);
    if(findArcSign(dec,arc).reversed){
        return findArcTail(dec,arc);
    }else{
        return findArcHead(dec,arc);
    }
}
static spqr_node findEffectiveArcTail(MATRECNetworkDecomposition *dec, spqr_arc arc){
    assert(dec);
    if(findArcSign(dec,arc).reversed){
        return findArcHead(dec,arc);
    }else{
        return findArcTail(dec,arc);
    }
}
static spqr_node findEffectiveArcHeadNoCompression(const MATRECNetworkDecomposition *dec, spqr_arc arc){
    assert(dec);
    if(findArcSignNoCompression(dec,arc).reversed){
        return findArcTailNoCompression(dec,arc);
    }else{
        return findArcHeadNoCompression(dec,arc);
    }
}
static spqr_node findEffectiveArcTailNoCompression(const MATRECNetworkDecomposition *dec, spqr_arc arc){
    assert(dec);
    if(findArcSignNoCompression(dec,arc).reversed){
        return findArcHeadNoCompression(dec,arc);
    }else{
        return findArcTailNoCompression(dec,arc);
    }
}
///Merge for signed union find of the arc directions.
///Is not symmetric, in the sense that the arc directions of coponent first are guaranteed not to change but those of second may change
///Based on whether one wants the reflection or not
static spqr_arc mergeArcSigns(MATRECNetworkDecomposition *dec, spqr_arc first, spqr_arc second,
                              bool reflectRelative) {
    assert(dec);
    assert(arcIsRepresentative(dec, first));
    assert(arcIsRepresentative(dec, second));
    assert(first != second); //We cannot merge a member into itself
    assert(first < dec->memArcs);
    assert(second < dec->memArcs);

    //The rank is stored as a negative number: we decrement it making the negative number larger.
    // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
    spqr_member firstRank = dec->arcs[first].representative;
    spqr_member secondRank = dec->arcs[second].representative;

    if (firstRank > secondRank) {
        swap_ints(&first, &second);
    }
    dec->arcs[second].representative = first;
    if (firstRank == secondRank) {
        --dec->arcs[first].representative;
    }
    //These boolean formula's cover all 16 possible cases, such that the relative orientation of the first is not changed
    bool equal = dec->arcs[first].reversed == dec->arcs[second].reversed;
    dec->arcs[second].reversed = (equal == reflectRelative);
    if(firstRank > secondRank){
        dec->arcs[first].reversed = (dec->arcs[first].reversed != reflectRelative);
    }
    return first;
}

static bool arcIsReversedNonRigid(const MATRECNetworkDecomposition *dec, spqr_arc arc){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    return dec->arcs[arc].reversed;
}


static spqr_element arcGetElement(const MATRECNetworkDecomposition * dec, spqr_arc arc){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);

    return dec->arcs[arc].element;
}
bool MATRECNetworkDecompositionContainsRow(const MATRECNetworkDecomposition * dec, MATREC_row row){
    assert(MATRECrowIsValid(row) && (int) row < dec->memRows);
    assert(dec);
    return SPQRarcIsValid(dec->rowArcs[row]);
}
bool MATRECNetworkDecompositionContainsColumn(const MATRECNetworkDecomposition *dec, MATREC_col col){
    assert(MATRECcolIsValid(col) && (int) col < dec->memColumns);
    assert(dec);
    return SPQRarcIsValid(dec->columnArcs[col]);
}
static void setDecompositionColumnArc(MATRECNetworkDecomposition *dec, MATREC_col col, spqr_arc arc){
    assert(MATRECcolIsValid(col) && (int)col < dec->memColumns);
    assert(dec);
    assert(SPQRarcIsValid(arc));
    dec->columnArcs[col] = arc;
}
static void setDecompositionRowArc(MATRECNetworkDecomposition *dec, MATREC_row row, spqr_arc arc){
    assert(MATRECrowIsValid(row) && (int) row < dec->memRows);
    assert(dec);
    assert(SPQRarcIsValid(arc));
    dec->rowArcs[row] = arc;
}
static spqr_arc getDecompositionColumnArc(const MATRECNetworkDecomposition *dec, MATREC_col col){
    assert(MATRECcolIsValid(col) && (int) col < dec->memColumns);
    assert(dec);
    return dec->columnArcs[col];
}
static spqr_arc getDecompositionRowArc(const MATRECNetworkDecomposition *dec, MATREC_row row){
    assert(MATRECrowIsValid(row) && (int) row < dec->memRows);
    assert(dec);
    return dec->rowArcs[row];
}

MATREC_ERROR MATRECNetworkDecompositionCreate(MATREC * env, MATRECNetworkDecomposition **pDecomposition, int numRows, int numColumns){
    assert(env);
    assert(pDecomposition);
    assert(!*pDecomposition);

    MATREC_CALL(MATRECallocBlock(env, pDecomposition));
    MATRECNetworkDecomposition *dec = *pDecomposition;
    dec->env = env;

    //Initialize arc array data
    int initialMemArcs = 8;
    {
        assert(initialMemArcs > 0);
        dec->memArcs = initialMemArcs;
        dec->numArcs = 0;
        MATREC_CALL(MATRECallocBlockArray(env, &dec->arcs, (size_t) dec->memArcs));
        for (spqr_arc i = 0; i < dec->memArcs; ++i) {
            dec->arcs[i].arcListNode.next = i + 1;
            dec->arcs[i].member = SPQR_INVALID_MEMBER;
        }
        dec->arcs[dec->memArcs - 1].arcListNode.next = SPQR_INVALID_ARC;
        dec->firstFreeArc = 0;
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
        MATREC_CALL(MATRECallocBlockArray(env, &dec->rowArcs, (size_t) dec->memRows));
        for (int i = 0; i < dec->memRows; ++i) {
            dec->rowArcs[i] = SPQR_INVALID_ARC;
        }
    }
    //Initialize mappings for columns
    {
        dec->memColumns = numColumns;
        dec->numColumns = 0;
        MATREC_CALL(MATRECallocBlockArray(env, &dec->columnArcs, (size_t) dec->memColumns));
        for (int i = 0; i < dec->memColumns; ++i) {
            dec->columnArcs[i] = SPQR_INVALID_ARC;
        }
    }

    dec->numConnectedComponents = 0;
    return MATREC_OKAY;
}

void MATRECNetworkDecompositionFree(MATRECNetworkDecomposition **pDec){
    assert(pDec);
    assert(*pDec);

    MATRECNetworkDecomposition *dec = *pDec;
    MATRECfreeBlockArray(dec->env, &dec->columnArcs);
    MATRECfreeBlockArray(dec->env, &dec->rowArcs);
    MATRECfreeBlockArray(dec->env, &dec->nodes);
    MATRECfreeBlockArray(dec->env, &dec->members);
    MATRECfreeBlockArray(dec->env, &dec->arcs);

    MATRECfreeBlock(dec->env, pDec);

}
static spqr_arc getFirstMemberArc(const MATRECNetworkDecomposition * dec, spqr_member member){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    return dec->members[member].firstArc;
}
static spqr_arc getNextMemberArc(const MATRECNetworkDecomposition * dec, spqr_arc arc){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    arc = dec->arcs[arc].arcListNode.next;
    return arc;
}
static spqr_arc getPreviousMemberArc(const MATRECNetworkDecomposition *dec, spqr_arc arc){
    assert(dec);
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    arc = dec->arcs[arc].arcListNode.previous;
    return arc;
}

static void addArcToMemberArcList(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_member member){
    spqr_arc firstMemberArc = getFirstMemberArc(dec, member);

    if(SPQRarcIsValid(firstMemberArc)){
        spqr_arc lastMemberArc = getPreviousMemberArc(dec, firstMemberArc);
        dec->arcs[arc].arcListNode.next = firstMemberArc;
        dec->arcs[arc].arcListNode.previous = lastMemberArc;
        dec->arcs[firstMemberArc].arcListNode.previous = arc;
        dec->arcs[lastMemberArc].arcListNode.next = arc;
    }else{
        assert(dec->members[member].numArcs == 0);
        dec->arcs[arc].arcListNode.next = arc;
        dec->arcs[arc].arcListNode.previous = arc;
    }
    dec->members[member].firstArc = arc;//TODO: update this in case of row/column arcs to make memory ordering nicer?
    ++(dec->members[member].numArcs);
}
static MATREC_ERROR createArc(MATRECNetworkDecomposition *dec, spqr_member member,bool reversed, spqr_arc *pArc) {
    assert(dec);
    assert(pArc);
    assert(SPQRmemberIsInvalid(member) || memberIsRepresentative(dec, member));

    spqr_arc index = dec->firstFreeArc;
    if (SPQRarcIsValid(index)) {
        dec->firstFreeArc = dec->arcs[index].arcListNode.next;
    } else {
        //Enlarge array, no free nodes in arc list
        int newSize = 2 * dec->memArcs;
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &dec->arcs, (size_t) newSize));
        for (int i = dec->memArcs + 1; i < newSize; ++i) {
            dec->arcs[i].arcListNode.next = i + 1;
            dec->arcs[i].member = SPQR_INVALID_MEMBER;
        }
        dec->arcs[newSize - 1].arcListNode.next = SPQR_INVALID_ARC;
        dec->firstFreeArc = dec->memArcs + 1;
        index = dec->memArcs;
        dec->memArcs = newSize;
    }
    //TODO: Is defaulting these here necessary?
    dec->arcs[index].tail = SPQR_INVALID_NODE;
    dec->arcs[index].head = SPQR_INVALID_NODE;
    dec->arcs[index].member = member;
    dec->arcs[index].childMember = SPQR_INVALID_MEMBER;
    dec->arcs[index].reversed = reversed;

    dec->arcs[index].headArcListNode.next = SPQR_INVALID_ARC;
    dec->arcs[index].headArcListNode.previous = SPQR_INVALID_ARC;
    dec->arcs[index].tailArcListNode.next = SPQR_INVALID_ARC;
    dec->arcs[index].tailArcListNode.previous = SPQR_INVALID_ARC;

    dec->numArcs++;

    *pArc = index;

    return MATREC_OKAY;
}
static MATREC_ERROR createRowArc(MATRECNetworkDecomposition *dec, spqr_member member, spqr_arc *pArc,
                                 MATREC_row row, bool reversed){
    MATREC_CALL(createArc(dec,member,reversed,pArc));
    setDecompositionRowArc(dec,row,*pArc);
    addArcToMemberArcList(dec,*pArc,member);
    dec->arcs[*pArc].element = MATRECrowToElement(row);

    return MATREC_OKAY;
}
static MATREC_ERROR createColumnArc(MATRECNetworkDecomposition *dec, spqr_member member, spqr_arc *pArc,
                                    MATREC_col column, bool reversed){
    MATREC_CALL(createArc(dec,member,reversed,pArc));
    setDecompositionColumnArc(dec,column,*pArc);
    addArcToMemberArcList(dec,*pArc,member);
    dec->arcs[*pArc].element = MATRECcolumnToElement(column);

    return MATREC_OKAY;
}
static MATREC_ERROR createMember(MATRECNetworkDecomposition *dec, SPQRMemberType type, spqr_member * pMember){
    assert(dec);
    assert(pMember);

    if(dec->numMembers == dec->memMembers){
        dec->memMembers *= 2;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&dec->members,(size_t) dec->memMembers));
    }
    MATRECNetworkDecompositionMember *data = &dec->members[dec->numMembers];
    data->markerOfParent = SPQR_INVALID_ARC;
    data->markerToParent = SPQR_INVALID_ARC;
    data->firstArc = SPQR_INVALID_ARC;
    data->representativeMember = SPQR_INVALID_MEMBER;
    data->numArcs = 0;
    data->parentMember = SPQR_INVALID_MEMBER;
    data->type = type;

    *pMember = dec->numMembers;

    dec->numMembers++;
    return MATREC_OKAY;
}

static MATREC_ERROR createNode(MATRECNetworkDecomposition *dec, spqr_node * pNode){

    if(dec->numNodes == dec->memNodes){
        dec->memNodes*=2;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&dec->nodes,(size_t) dec->memNodes));
    }
    *pNode = dec->numNodes;
    dec->nodes[dec->numNodes].representativeNode = SPQR_INVALID_NODE;
    dec->nodes[dec->numNodes].firstArc = SPQR_INVALID_ARC;
    dec->nodes[dec->numNodes].numArcs = 0;
    dec->numNodes++;

    return MATREC_OKAY;
}
static void removeArcFromNodeArcList(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_node node, bool nodeIsHead){
    MATRECNetworkDecompositionArcListNode * arcListNode = nodeIsHead ? &dec->arcs[arc].headArcListNode : &dec->arcs[arc].tailArcListNode;

    if(dec->nodes[node].numArcs == 1){
        dec->nodes[node].firstArc = SPQR_INVALID_ARC;
    }else{
        spqr_arc next_arc = arcListNode->next;
        spqr_arc prev_arc = arcListNode->previous;
        MATRECNetworkDecompositionArcListNode * nextListNode = findArcHead(dec, next_arc) == node ? &dec->arcs[next_arc].headArcListNode : &dec->arcs[next_arc].tailArcListNode;//TODO: finds necessary?
        MATRECNetworkDecompositionArcListNode * prevListNode = findArcHead(dec, prev_arc) == node ? &dec->arcs[prev_arc].headArcListNode : &dec->arcs[prev_arc].tailArcListNode;//TODO: finds necessary?

        nextListNode->previous = prev_arc;
        prevListNode->next = next_arc;

        if(dec->nodes[node].firstArc == arc){
            dec->nodes[node].firstArc = next_arc; //TODO: fix this if we want fixed ordering for tree/nontree arcs in memory
        }
    }
    //TODO: empty arcListNode's data? Might not be all that relevant
    --(dec->nodes[node].numArcs);
}
static void addArcToNodeArcList(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_node node, bool nodeIsHead){
    assert(nodeIsRepresentative(dec,node));

    spqr_arc firstNodeArc = getFirstNodeArc(dec, node);

    MATRECNetworkDecompositionArcListNode * arcListNode = nodeIsHead ? &dec->arcs[arc].headArcListNode : &dec->arcs[arc].tailArcListNode;
    if(SPQRarcIsValid(firstNodeArc)){
        bool nextIsHead = findArcHead(dec,firstNodeArc) == node;
        MATRECNetworkDecompositionArcListNode *nextListNode = nextIsHead ? &dec->arcs[firstNodeArc].headArcListNode : &dec->arcs[firstNodeArc].tailArcListNode;
        spqr_arc lastNodeArc = nextListNode->previous;

        arcListNode->next = firstNodeArc;
        arcListNode->previous = lastNodeArc;


        bool previousIsHead = findArcHead(dec,lastNodeArc) == node;
        MATRECNetworkDecompositionArcListNode *previousListNode = previousIsHead ? &dec->arcs[lastNodeArc].headArcListNode : &dec->arcs[lastNodeArc].tailArcListNode;
        previousListNode->next = arc;
        nextListNode->previous = arc;

    }else{
        arcListNode->next = arc;
        arcListNode->previous = arc;
    }
    dec->nodes[node].firstArc = arc; //TODO: update this in case of row/column arcs to make memory ordering nicer?er?
    ++dec->nodes[node].numArcs;
    if(nodeIsHead){
        dec->arcs[arc].head = node;
    }else{
        dec->arcs[arc].tail = node;
    }
}
static void setArcHeadAndTail(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_node head, spqr_node tail){
    addArcToNodeArcList(dec,arc,head,true);
    addArcToNodeArcList(dec,arc,tail,false);
}
static void clearArcHeadAndTail(MATRECNetworkDecomposition *dec, spqr_arc arc){
    removeArcFromNodeArcList(dec,arc,findArcHead(dec,arc),true);
    removeArcFromNodeArcList(dec,arc,findArcTail(dec,arc),false);
    dec->arcs[arc].head = SPQR_INVALID_NODE;
    dec->arcs[arc].tail = SPQR_INVALID_NODE;
}
static void changeArcHead(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_node oldHead, spqr_node newHead){
    assert(nodeIsRepresentative(dec,oldHead));
    assert(nodeIsRepresentative(dec,newHead));
    removeArcFromNodeArcList(dec,arc,oldHead,true);
    addArcToNodeArcList(dec,arc,newHead,true);
}
static void changeArcTail(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_node oldTail, spqr_node newTail){
    assert(nodeIsRepresentative(dec,oldTail));
    assert(nodeIsRepresentative(dec,newTail));
    removeArcFromNodeArcList(dec,arc,oldTail,false);
    addArcToNodeArcList(dec,arc,newTail,false);
}

static int nodeDegree(MATRECNetworkDecomposition *dec, spqr_node node){
    assert(dec);
    assert(SPQRnodeIsValid(node));
    assert(node < dec->memNodes);
    return dec->nodes[node].numArcs;
}
static SPQRMemberType getMemberType(const MATRECNetworkDecomposition *dec, spqr_member member){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));
    return dec->members[member].type;
}
static void updateMemberType(const MATRECNetworkDecomposition *dec, spqr_member member, SPQRMemberType type){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));

    dec->members[member].type = type;
}
static spqr_arc markerToParent(const MATRECNetworkDecomposition *dec, spqr_member member){
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

static void updateMemberParentInformation(MATRECNetworkDecomposition *dec, const spqr_member newMember, const spqr_member toRemove){
    assert(memberIsRepresentative(dec,newMember));
    assert(findMemberNoCompression(dec,toRemove) == newMember);

    dec->members[newMember].markerOfParent = dec->members[toRemove].markerOfParent;
    dec->members[newMember].markerToParent = dec->members[toRemove].markerToParent;
    dec->members[newMember].parentMember = dec->members[toRemove].parentMember;

    dec->members[toRemove].markerOfParent = SPQR_INVALID_ARC;
    dec->members[toRemove].markerToParent = SPQR_INVALID_ARC;
    dec->members[toRemove].parentMember = SPQR_INVALID_MEMBER;
}
static spqr_arc markerOfParent(const MATRECNetworkDecomposition *dec, spqr_member member) {
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));
    return dec->members[member].markerOfParent;
}



static int getNumMemberArcs(const MATRECNetworkDecomposition * dec, spqr_member member){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(memberIsRepresentative(dec,member));
    return dec->members[member].numArcs;
}

static int getNumNodes(const MATRECNetworkDecomposition *dec){
    assert(dec);
    return dec->numNodes;
}
static int getNumMembers(const MATRECNetworkDecomposition *dec){
    assert(dec);
    return dec->numMembers;
}
static MATREC_ERROR createStandaloneParallel(MATRECNetworkDecomposition *dec, MATREC_col * columns, bool * reversed,
                                             int num_columns, MATREC_row row, spqr_member * pMember){
    spqr_member member;
    MATREC_CALL(createMember(dec, num_columns <= 1 ? SPQR_MEMBERTYPE_LOOP : SPQR_MEMBERTYPE_PARALLEL, &member));

    spqr_arc row_arc;
    MATREC_CALL(createRowArc(dec,member,&row_arc,row,num_columns <= 1));

    spqr_arc col_arc;
    for (int i = 0; i < num_columns; ++i) {
        MATREC_CALL(createColumnArc(dec,member,&col_arc,columns[i],reversed[i]));
    }
    *pMember = member;

    ++dec->numConnectedComponents;
    //TODO: double check that this function is not used to create connected parallel (We could also compute this quantity by tracking # of representative members - number of marker pairs
    return MATREC_OKAY;
}

//TODO: fix tracking connectivity more cleanly, should not be left up to the algorithms ideally
static MATREC_ERROR createConnectedParallel(MATRECNetworkDecomposition *dec, MATREC_col * columns, bool * reversed,
                                            int num_columns, MATREC_row row, spqr_member * pMember){
    spqr_member member;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &member));

    spqr_arc row_arc;
    MATREC_CALL(createRowArc(dec,member,&row_arc,row,false));

    spqr_arc col_arc;
    for (int i = 0; i < num_columns; ++i) {
        MATREC_CALL(createColumnArc(dec,member,&col_arc,columns[i],reversed[i]));
    }
    *pMember = member;

    return MATREC_OKAY;
}

static MATREC_ERROR createStandaloneSeries(MATRECNetworkDecomposition *dec, MATREC_row * rows, bool *reversed,
                                           int numRows, MATREC_col col, spqr_member * pMember){
    spqr_member member;
    MATREC_CALL(createMember(dec, numRows <= 1 ? SPQR_MEMBERTYPE_LOOP : SPQR_MEMBERTYPE_SERIES, &member));

    spqr_arc colArc;
    MATREC_CALL(createColumnArc(dec,member,&colArc,col,false));

    spqr_arc rowArc;
    for (int i = 0; i < numRows; ++i) {
        MATREC_CALL(createRowArc(dec,member,&rowArc,rows[i],!reversed[i]));
    }
    *pMember = member;
    ++dec->numConnectedComponents;
    return MATREC_OKAY;
}
static MATREC_ERROR createConnectedSeries(MATRECNetworkDecomposition *dec, MATREC_row * rows, bool *reversed,
                                          int numRows, MATREC_col col, spqr_member * pMember){
    spqr_member member;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &member));//TODO: check type here if numRows <= 1?

    spqr_arc colArc;
    MATREC_CALL(createColumnArc(dec,member,&colArc,col,false));

    spqr_arc rowArc;
    for (int i = 0; i < numRows; ++i) {
        MATREC_CALL(createRowArc(dec,member,&rowArc,rows[i],!reversed[i]));
    }
    *pMember = member;
    return MATREC_OKAY;
}

static void removeArcFromMemberArcList(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_member member){
    assert(findArcMemberNoCompression(dec,arc) == member);
    assert(memberIsRepresentative(dec,member));

    if(dec->members[member].numArcs == 1){
        dec->members[member].firstArc = SPQR_INVALID_ARC;

        //TODO: also set arcListNode to invalid, maybe? Not necessary probably
    }else{
        spqr_arc nextArc = dec->arcs[arc].arcListNode.next;
        spqr_arc prevArc = dec->arcs[arc].arcListNode.previous;

        dec->arcs[nextArc].arcListNode.previous = prevArc;
        dec->arcs[prevArc].arcListNode.next = nextArc;

        if(dec->members[member].firstArc == arc){
            dec->members[member].firstArc = nextArc; //TODO: fix this if we want fixed ordering for tree/nontree arcs in memory
        }
    }


    --(dec->members[member].numArcs);
}

typedef struct {
    spqr_arc arc;
    bool reversed;
} FindCycleCall;
static void process_arc(MATREC_row * fundamental_cycle_arcs, int * num_cycle_arcs,
                        FindCycleCall * callStack,
                        int * callStackSize,
                        spqr_arc arc,
                        const MATRECNetworkDecomposition * dec,
                        bool * fundamental_cycle_direction,
                        bool arcIsReversed){
    assert(arcIsTree(dec,arc));
    if(!arcIsMarker(dec,arc)){
        spqr_member current_member = findArcMemberNoCompression(dec, arc);
        if(markerToParent(dec,current_member) == arc){
            spqr_arc other_arc = markerOfParent(dec, current_member);
            assert(!arcIsTree(dec,other_arc));
            callStack[*callStackSize].arc = other_arc;
            callStack[*callStackSize].reversed = arcIsReversed;
            ++(*callStackSize);
        }else{
            spqr_element element = arcGetElement(dec,arc);
            assert(SPQRelementIsRow(element));
            MATREC_row row = SPQRelementToRow(element);
            fundamental_cycle_arcs[*num_cycle_arcs] = row;
            fundamental_cycle_direction[*num_cycle_arcs] = arcIsReversed;
            ++(*num_cycle_arcs);
        }
    }else{
        spqr_member child_member = findArcChildMemberNoCompression(dec, arc);
        spqr_arc other_arc = markerToParent(dec, child_member);
        assert(!arcIsTree(dec,other_arc));
        callStack[*callStackSize].arc = other_arc;
        callStack[*callStackSize].reversed = arcIsReversed;
        ++(*callStackSize);
    }
}

static int decompositionGetFundamentalCycleRows(const MATRECNetworkDecomposition *dec, MATREC_col column, MATREC_row * output,
                                                bool * computedSignStorage){
    spqr_arc arc = getDecompositionColumnArc(dec, column);
    if(SPQRarcIsInvalid(arc)){
        return 0;
    }
    int num_rows = 0;

    FindCycleCall * callStack;
    //TODO: probably an overkill amount of memory allocated here... How can we allocate just enough?
    MATREC_ERROR result = MATRECallocBlockArray(dec->env,&callStack,(size_t) dec->memRows);
    if(result != MATREC_OKAY){
        return -1;
    }
    int callStackSize = 1;
    callStack[0].arc = arc;
    callStack[0].reversed = false; //TODO: check?

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
        spqr_arc nodeArc;
    } DFSCallData;
    DFSCallData * pathSearchCallStack;
    result = MATRECallocBlockArray(dec->env,&pathSearchCallStack,(size_t) dec->numNodes);
    if(result != MATREC_OKAY){
        return -1;
    }
    int pathSearchCallStackSize = 0;

    while(callStackSize > 0){
        spqr_arc column_arc = callStack[callStackSize - 1].arc;
        bool reverseEverything = callStack[callStackSize-1].reversed;
        --callStackSize;
        spqr_member column_arc_member = findArcMemberNoCompression(dec, column_arc);
        switch(getMemberType(dec,column_arc_member)){
            case SPQR_MEMBERTYPE_RIGID:
            {
                spqr_node source = findEffectiveArcTailNoCompression(dec, column_arc);
                spqr_node target = findEffectiveArcHeadNoCompression(dec, column_arc);

                assert(pathSearchCallStackSize == 0);
                pathSearchCallStack[0].node = source;
                pathSearchCallStack[0].nodeArc = getFirstNodeArc(dec,source);
                pathSearchCallStackSize++;
                while(pathSearchCallStackSize > 0){
                    DFSCallData * dfsData  = &pathSearchCallStack[pathSearchCallStackSize-1];
                    nodeVisited[dfsData->node] = true;
                    //cannot be a tree arc which is its parent
                    if(arcIsTree(dec,dfsData->nodeArc) &&
                       (pathSearchCallStackSize <= 1 || dfsData->nodeArc != pathSearchCallStack[pathSearchCallStackSize-2].nodeArc)){
                        spqr_node head = findEffectiveArcHeadNoCompression(dec, dfsData->nodeArc);
                        spqr_node tail = findEffectiveArcTailNoCompression(dec, dfsData->nodeArc);
                        spqr_node other = head == dfsData->node ? tail : head;
                        assert(other != dfsData->node);
                        assert(!nodeVisited[other]);
                        if(other == target){
                            break;
                        }
                        //We go up a level: add new node to the call stack

                        pathSearchCallStack[pathSearchCallStackSize].node = other;
                        pathSearchCallStack[pathSearchCallStackSize].nodeArc = getFirstNodeArc(dec,other);
                        ++pathSearchCallStackSize;
                        continue;
                    }
                    do{
                        dfsData->nodeArc = getNextNodeArcNoCompression(dec,dfsData->nodeArc,dfsData->node);
                        if(dfsData->nodeArc == getFirstNodeArc(dec,dfsData->node)){
                            --pathSearchCallStackSize;
                            dfsData = &pathSearchCallStack[pathSearchCallStackSize-1];
                        }else{
                            break;
                        }
                    }while(pathSearchCallStackSize > 0);
                }
                for (int i = 0; i < pathSearchCallStackSize; ++i) {
                    if(arcIsTree(dec,pathSearchCallStack[i].nodeArc)){
                        bool arcReversedInPath = findEffectiveArcHeadNoCompression(dec,pathSearchCallStack[i].nodeArc) == pathSearchCallStack[i].node;
                        //TODO: also check 'reversed'
                        process_arc(output,&num_rows,callStack,&callStackSize,pathSearchCallStack[i].nodeArc,dec,
                                    computedSignStorage,arcReversedInPath != reverseEverything);
                    }
                }

                pathSearchCallStackSize = 0;
                break;
            }
            case SPQR_MEMBERTYPE_PARALLEL:
            {
                bool columnReversed = arcIsReversedNonRigid(dec,column_arc);

                spqr_arc first_arc = getFirstMemberArc(dec, column_arc_member);
                spqr_arc iter_arc = first_arc;
                int tree_count = 0;
                do
                {
                    if(arcIsTree(dec,iter_arc)){
                        bool treeIsReversed = arcIsReversedNonRigid(dec,iter_arc);
                        process_arc(output,&num_rows,callStack,&callStackSize,iter_arc,dec,
                                    computedSignStorage,(columnReversed != treeIsReversed) != reverseEverything);
                        tree_count++;
                    }
                    iter_arc = getNextMemberArc(dec,iter_arc);
                }
                while(iter_arc != first_arc);
                if(tree_count != 1){
                    return -1;
                }
                break;
            }
            case SPQR_MEMBERTYPE_LOOP:
            case SPQR_MEMBERTYPE_SERIES:
            {
                bool columnReversed = arcIsReversedNonRigid(dec,column_arc);
                spqr_arc first_arc = getFirstMemberArc(dec, column_arc_member);
                spqr_arc iter_arc = first_arc;
                int nontree_count = 0;
                do
                {
                    if(arcIsTree(dec,iter_arc)){
                        bool treeIsReversed = arcIsReversedNonRigid(dec,iter_arc);
                        process_arc(output,&num_rows,callStack,&callStackSize,iter_arc,dec,
                                    computedSignStorage,(columnReversed == treeIsReversed) != reverseEverything);
                    }else{
                        nontree_count++;
                    }
                    iter_arc = getNextMemberArc(dec,iter_arc);
                }
                while(iter_arc != first_arc);
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
typedef struct{
    MATREC_row row;
    bool reversed;
} Nonzero;
static int qsort_comparison (const void * a, const void * b)
{
    Nonzero *s1 = (Nonzero *)a;
    Nonzero *s2 = (Nonzero *)b;

    if(s1->row > s2->row) {
        return 1;
    }
    else if(s1->row == s2->row) {
        return 0;
    }
    else {
        return -1;
    }
}
bool MATRECNetworkDecompositionVerifyCycle(const MATRECNetworkDecomposition * dec, MATREC_col column, MATREC_row * column_rows,
                                           double * column_values, int num_rows,
                                           MATREC_row * computed_column_storage,
                                           bool * computedSignStorage){
    int num_found_rows = decompositionGetFundamentalCycleRows(dec,column,computed_column_storage,computedSignStorage);

    if(num_found_rows != num_rows){
        return false;
    }
    if(num_rows == 0){
        return true;
    }

    Nonzero array[num_rows];
    for (int i = 0; i < num_rows; ++i) {
        array[i].row = computed_column_storage[i];
        array[i].reversed = computedSignStorage[i];
    }
    qsort(array,(size_t) num_rows,sizeof(Nonzero),qsort_comparison);

    Nonzero secondArray[num_rows];
    for (int i = 0; i < num_rows; ++i) {
        secondArray[i].row = column_rows[i];
        secondArray[i].reversed = column_values[i] < 0.0;
    }

    qsort(secondArray,(size_t) num_rows,sizeof(Nonzero),qsort_comparison);

    for (int i = 0; i < num_rows; ++i) {
        if(array[i].row != secondArray[i].row){
            return false;
        }
        if(array[i].reversed != secondArray[i].reversed){
            return false;
        }
    }
    return true;
}

static spqr_member largestMemberID(const MATRECNetworkDecomposition *dec){
    return dec->numMembers;
}
static spqr_arc largestArcID(const MATRECNetworkDecomposition *dec){
    return dec->numArcs;
}
static spqr_node largestNodeID(const MATRECNetworkDecomposition *dec){
    return dec->numNodes;
}
static int numConnectedComponents(const MATRECNetworkDecomposition *dec){
    return dec->numConnectedComponents;
}
static MATREC_ERROR createChildMarker(MATRECNetworkDecomposition *dec, spqr_member member, spqr_member child, bool isTree,
                                    spqr_arc * pArc, bool reversed){
    MATREC_CALL(createArc(dec,member,reversed,pArc)); //TODO: fix
    dec->arcs[*pArc].element = isTree ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;
    dec->arcs[*pArc].childMember = child;

    addArcToMemberArcList(dec,*pArc,member);
    return MATREC_OKAY;
}
static MATREC_ERROR createParentMarker(MATRECNetworkDecomposition *dec, spqr_member member, bool isTree, spqr_member parent, spqr_arc parentMarker
        , spqr_arc * arc,bool reversed){

    MATREC_CALL(createArc(dec,member,reversed,arc)); //TODO: fix
    dec->arcs[*arc].element = isTree ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;

    addArcToMemberArcList(dec,*arc,member);

    dec->members[member].parentMember = parent;
    dec->members[member].markerOfParent = parentMarker;
    dec->members[member].markerToParent = *arc;
    return MATREC_OKAY;
}
static MATREC_ERROR createMarkerPair(MATRECNetworkDecomposition *dec, spqr_member parentMember, spqr_member childMember,
                                   bool parentIsTree,
                                   bool childMarkerReversed,
                                   bool parentMarkerReversed){
    spqr_arc parentToChildMarker = SPQR_INVALID_ARC;
    MATREC_CALL(createChildMarker(dec,parentMember,childMember,parentIsTree,&parentToChildMarker,childMarkerReversed));

    spqr_arc childToParentMarker = SPQR_INVALID_ARC;
    MATREC_CALL(createParentMarker(dec,childMember,!parentIsTree,parentMember,parentToChildMarker,&childToParentMarker,parentMarkerReversed));

    return MATREC_OKAY;
}
static MATREC_ERROR createMarkerPairWithReferences(MATRECNetworkDecomposition *dec, spqr_member parentMember, spqr_member childMember, bool parentIsTree,
                                                 bool childMarkerReversed,
                                                 bool parentMarkerReversed,
                                                 spqr_arc * parentToChild, spqr_arc *childToParent){
    MATREC_CALL(createChildMarker(dec,parentMember,childMember,parentIsTree,parentToChild,childMarkerReversed));
    MATREC_CALL(createParentMarker(dec,childMember,!parentIsTree,parentMember,*parentToChild,childToParent,parentMarkerReversed));

    return MATREC_OKAY;
}

static void moveArcToNewMember(MATRECNetworkDecomposition *dec, spqr_arc arc, spqr_member oldMember, spqr_member newMember){
    assert(SPQRarcIsValid(arc));
    assert(arc < dec->memArcs);
    assert(dec);

    assert(memberIsRepresentative(dec,oldMember));
    assert(memberIsRepresentative(dec,newMember));
    //Need to change the arc's member, remove it from the current member list and add it to the new member list
    assert(findArcMemberNoCompression(dec,arc) == oldMember);

    removeArcFromMemberArcList(dec,arc,oldMember);
    addArcToMemberArcList(dec,arc,newMember);

    dec->arcs[arc].member = newMember;

    //If this arc has a childMember, update the information correctly!
    spqr_member childMember = dec->arcs[arc].childMember;
    if(SPQRmemberIsValid(childMember)){
        spqr_member childRepresentative = findArcChildMember(dec, arc);
        dec->members[childRepresentative].parentMember = newMember;
    }
    //If this arc is a marker to the parent, update the child arc marker of the parent to reflect the move
    if(dec->members[oldMember].markerToParent == arc){
        dec->members[newMember].markerToParent = arc;
        dec->members[newMember].parentMember = dec->members[oldMember].parentMember;
        dec->members[newMember].markerOfParent = dec->members[oldMember].markerOfParent;

        assert(findArcChildMemberNoCompression(dec,dec->members[oldMember].markerOfParent) == oldMember);
        dec->arcs[dec->members[oldMember].markerOfParent].childMember = newMember;
    }
}
static void mergeMemberArcList(MATRECNetworkDecomposition *dec, spqr_member toMergeInto, spqr_member toRemove){
    spqr_arc firstIntoArc = getFirstMemberArc(dec, toMergeInto);
    spqr_arc firstFromArc = getFirstMemberArc(dec, toRemove);
    assert(SPQRarcIsValid(firstIntoArc));
    assert(SPQRarcIsValid(firstFromArc));

    spqr_arc lastIntoArc = getPreviousMemberArc(dec, firstIntoArc);
    spqr_arc lastFromArc = getPreviousMemberArc(dec, firstFromArc);

    //Relink linked lists to merge them effectively
    dec->arcs[firstIntoArc].arcListNode.previous = lastFromArc;
    dec->arcs[lastIntoArc].arcListNode.next = firstFromArc;
    dec->arcs[firstFromArc].arcListNode.previous = lastIntoArc;
    dec->arcs[lastFromArc].arcListNode.next = firstIntoArc;

    //Clean up old
    dec->members[toMergeInto].numArcs += dec->members[toRemove].numArcs;
    dec->members[toRemove].numArcs = 0;
    dec->members[toRemove].firstArc = SPQR_INVALID_ARC;

}

static void changeLoopToSeries(MATRECNetworkDecomposition * dec, spqr_member member){
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(dec);
    assert((getMemberType(dec,member) == SPQR_MEMBERTYPE_PARALLEL || getMemberType(dec, member) == SPQR_MEMBERTYPE_SERIES ||
            getMemberType(dec,member) == SPQR_MEMBERTYPE_LOOP) && getNumMemberArcs(dec, member) == 2);
    assert(memberIsRepresentative(dec,member));
    dec->members[member].type = SPQR_MEMBERTYPE_SERIES;
}
static void changeLoopToParallel(MATRECNetworkDecomposition * dec, spqr_member member){
    assert(SPQRmemberIsValid(member));
    assert(member < dec->memMembers);
    assert(dec);
    assert((getMemberType(dec,member) == SPQR_MEMBERTYPE_PARALLEL || getMemberType(dec, member) == SPQR_MEMBERTYPE_SERIES ||
            getMemberType(dec,member) == SPQR_MEMBERTYPE_LOOP) && getNumMemberArcs(dec, member) == 2);
    assert(memberIsRepresentative(dec,member));
    dec->members[member].type = SPQR_MEMBERTYPE_PARALLEL;
}
bool MATRECNetworkDecompositionIsMinimal(const MATRECNetworkDecomposition * dec){
    //Relies on parents/children etc. being set correctly in the tree
    bool isMinimal = true;
    for (spqr_member member = 0; member < dec->numMembers; ++member) {
        if (!memberIsRepresentative(dec, member) || getMemberType(dec,member) == SPQR_MEMBERTYPE_UNASSIGNED ){
            continue;
        } //TODO: fix loop making INVALID members here... this is not a pretty way to do this
        spqr_member memberParent = findMemberParentNoCompression(dec, member);
        if(SPQRmemberIsValid(memberParent)){
            assert(memberIsRepresentative(dec,memberParent));
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

static void decreaseNumConnectedComponents(MATRECNetworkDecomposition *dec, int by){
    dec->numConnectedComponents-= by;
    assert(dec->numConnectedComponents >= 1);
}

static void reorderComponent(MATRECNetworkDecomposition *dec, spqr_member newRoot){
    assert(dec);
    assert(memberIsRepresentative(dec,newRoot));
    //If the newRoot has no parent, it is already the root, so then there's no need to reorder.
    if(SPQRmemberIsValid(dec->members[newRoot].parentMember)){
        spqr_member member = findMemberParent(dec, newRoot);
        spqr_member newParent = newRoot;
        spqr_arc newMarkerToParent = dec->members[newRoot].markerOfParent;
        spqr_arc markerOfNewParent = dec->members[newRoot].markerToParent;

        //Recursively update the parent
        do{
            assert(SPQRmemberIsValid(member));
            assert(SPQRmemberIsValid(newParent));
            spqr_member oldParent = findMemberParent(dec, member);
            spqr_arc oldMarkerToParent = dec->members[member].markerToParent;
            spqr_arc oldMarkerOfParent = dec->members[member].markerOfParent;

            dec->members[member].markerToParent = newMarkerToParent;
            dec->members[member].markerOfParent = markerOfNewParent;
            dec->members[member].parentMember = newParent;
            dec->arcs[markerOfNewParent].childMember = member;
            dec->arcs[newMarkerToParent].childMember = SPQR_INVALID_MEMBER;

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
        dec->members[newRoot].markerToParent = SPQR_INVALID_ARC;
        dec->members[newRoot].markerOfParent = SPQR_INVALID_ARC;

    }
}

static void arcToDot(FILE * stream, const MATRECNetworkDecomposition * dec,
                      spqr_arc arc, unsigned long dot_head, unsigned long dot_tail, bool useElementNames){
    assert(SPQRarcIsValid(arc));
    spqr_member member = findArcMemberNoCompression(dec, arc);
    SPQRMemberType member_type = getMemberType(dec, member);
    char type = typeToChar(member_type);
    const char* color = arcIsTree(dec,arc) ? ",color=red" :",color=blue";

    int arc_name = arc;

    if(markerToParent(dec,member) == arc){
        if(useElementNames){
            arc_name = -1;
        }
        fprintf(stream, "    %c_%d_%lu -> %c_p_%d [label=\"%d\",style=dashed%s];\n", type, member, dot_tail, type, member, arc_name, color);
        fprintf(stream, "    %c_p_%d -> %c_%d_%lu [label=\"%d\",style=dashed%s];\n", type, member, type, member, dot_head, arc_name, color);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
        fprintf(stream, "    %c_p_%d [style=dashed];\n", type, member);
    }else if(arcIsMarker(dec,arc)){
        spqr_member child = findArcChildMemberNoCompression(dec, arc);
        char childType = typeToChar(getMemberType(dec,child));
        if(useElementNames){
            arc_name = -1;
        }
        fprintf(stream, "    %c_%d_%lu -> %c_c_%d [label=\"%d\",style=dotted%s];\n", type, member, dot_tail, type, child, arc_name, color);
        fprintf(stream, "    %c_c_%d -> %c_%d_%lu [label=\"%d\",style=dotted%s];\n", type, child, type, member, dot_head, arc_name, color);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
        fprintf(stream, "    %c_c_%d [style=dotted];\n", type, child);
        fprintf(stream, "    %c_p_%d -> %c_c_%d [style=dashed,dir=forward];\n", childType, child, type, child);
    }else{
        if(useElementNames){
            spqr_element element = dec->arcs[arc].element;
            if(SPQRelementIsRow(element)){
                arc_name = (int) SPQRelementToRow(element);
            }else{
                arc_name = (int) SPQRelementToColumn(element);
            }
        }

        fprintf(stream, "    %c_%d_%lu -> %c_%d_%lu [label=\"%d \",style=bold%s];\n", type, member, dot_tail, type, member, dot_head,
                arc_name, color);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
        fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
    }
}

static void decompositionToDot(FILE * stream, const MATRECNetworkDecomposition *dec, bool useElementNames ){
    fprintf(stream, "//decomposition\ndigraph decomposition{\n   compound = true;\n");
    for (spqr_member member = 0; member < dec->numMembers; ++member){
        if(!memberIsRepresentative(dec,member)) continue;
        fprintf(stream,"   subgraph member_%d{\n",member);
        switch(getMemberType(dec,member)){
            case SPQR_MEMBERTYPE_RIGID:
            {
                spqr_arc first_arc = getFirstMemberArc(dec, member);
                spqr_arc arc = first_arc;
                do
                {
                    unsigned long arcHead = (unsigned long) findEffectiveArcHeadNoCompression(dec,arc);
                    unsigned long arcTail = (unsigned long) findEffectiveArcTailNoCompression(dec,arc);
                    arcToDot(stream,dec,arc,arcHead,arcTail,useElementNames);
                    arc = getNextMemberArc(dec,arc);
                }
                while(arc != first_arc);
                break;
            }
            case SPQR_MEMBERTYPE_LOOP:
            case SPQR_MEMBERTYPE_PARALLEL:
            {
                spqr_arc first_arc = getFirstMemberArc(dec, member);
                spqr_arc arc = first_arc;
                do
                {
                    if(arcIsReversedNonRigid(dec,arc)){
                        arcToDot(stream,dec,arc,1,0,useElementNames);
                    }else{
                        arcToDot(stream,dec,arc,0,1,useElementNames);
                    }
                    arc = getNextMemberArc(dec,arc);
                }
                while(arc != first_arc);
                break;
            }
            case SPQR_MEMBERTYPE_SERIES:
            {
                unsigned long i = 0;
                unsigned long num_member_arcs = (unsigned long) getNumMemberArcs(dec, member);
                spqr_arc first_arc = getFirstMemberArc(dec, member);
                spqr_arc arc = first_arc;
                do {
                    unsigned long head = i;
                    unsigned long tail = (i+1) % num_member_arcs;
                    if(arcIsReversedNonRigid(dec,arc)){
                        unsigned long temp = head;
                        head = tail;
                        tail = temp;
                    }
                    arcToDot(stream, dec, arc, head,tail,useElementNames);
                    arc = getNextMemberArc(dec, arc);
                    i++;
                } while (arc != first_arc);
                break;
            }
            case SPQR_MEMBERTYPE_UNASSIGNED:
                break;
        }
        fprintf(stream,"   }\n");
    }
    fprintf(stream,"}\n");
}

static MATREC_ERROR mergeGivenMemberIntoParent(MATRECNetworkDecomposition *dec,
                                        spqr_member member,
                                        spqr_member parent,
                                        spqr_arc parentToChild,
                                        spqr_arc childToParent,
                                        bool headToHead,
                                        spqr_member * mergedMember){
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));
    assert(SPQRmemberIsValid(parent));
    assert(memberIsRepresentative(dec,parent));
    assert(findMemberParentNoCompression(dec,member) == parent);
    assert(markerOfParent(dec, member) == parentToChild);
    assert(markerToParent(dec, member) == childToParent);

    removeArcFromMemberArcList(dec,parentToChild,parent);
    removeArcFromMemberArcList(dec,childToParent,member);

    //TODO: need to account for reversals here
    spqr_node parentArcNodes[2] = {findEffectiveArcTail(dec, parentToChild), findEffectiveArcHead(dec, parentToChild)};
    spqr_node childArcNodes[2] = {findEffectiveArcTail(dec, childToParent), findEffectiveArcHead(dec, childToParent)};

    clearArcHeadAndTail(dec,parentToChild);
    clearArcHeadAndTail(dec,childToParent);

    spqr_node first = childArcNodes[headToHead ? 0 : 1];
    spqr_node second = childArcNodes[headToHead ? 1 : 0];
    {
        spqr_node newNode = mergeNodes(dec, parentArcNodes[0], first);
        spqr_node toRemoveFrom = newNode == first ? parentArcNodes[0] : first;
        mergeNodeArcList(dec,newNode,toRemoveFrom);
    }
    {
        spqr_node newNode = mergeNodes(dec, parentArcNodes[1], second);
        spqr_node toRemoveFrom = newNode == second ? parentArcNodes[1] : second;
        mergeNodeArcList(dec,newNode,toRemoveFrom);
    }


    spqr_member newMember = mergeMembers(dec, member, parent);
    spqr_member toRemoveFrom = newMember == member ? parent : member;
    mergeMemberArcList(dec,newMember,toRemoveFrom);
    if(toRemoveFrom == parent){
        updateMemberParentInformation(dec,newMember,toRemoveFrom);
    }
    updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);
    *mergedMember = newMember;
    return MATREC_OKAY;
}

static int max(int a, int b){
    return (a > b) ? a : b;
}

typedef int path_arc_id;
#define INVALID_PATH_ARC (-1)

static bool pathArcIsInvalid(const path_arc_id arc) {
    return arc < 0;
}

static bool pathArcIsValid(const path_arc_id arc) {
    return !pathArcIsInvalid(arc);
}

typedef struct {
    spqr_arc arc;
    spqr_node arcHead; //These can be used in various places to prevent additional find()'s
    spqr_node arcTail;
    path_arc_id nextMember;
    path_arc_id nextOverall;
    bool reversed;
} PathArcListNode;

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
    REDUCEDMEMBER_TYPE_UNASSIGNED = 0,
    REDUCEDMEMBER_TYPE_CYCLE = 1,
    REDUCEDMEMBER_TYPE_MERGED = 2,
    REDUCEDMEMBER_TYPE_NOT_NETWORK = 3
    //TODO fix
} ReducedMemberType;

typedef enum {
    INTO_HEAD = 0,
    INTO_TAIL = 1,
    OUT_HEAD = 2,
    OUT_TAIL = 3
} MemberPathType;

static bool isInto(MemberPathType type){
    return type < 2;
}
static bool isHead(MemberPathType type){
    return (type & 1) == 0;
}

typedef struct {
    spqr_member member;
    spqr_member rootMember;
    int depth;
    ReducedMemberType type;
    reduced_member_id parent;

    children_idx firstChild;
    children_idx numChildren;

    path_arc_id firstPathArc;
    int numPathArcs;

    bool reverseArcs;
    spqr_node rigidPathStart;
    spqr_node rigidPathEnd;

    bool pathBackwards;

    int numPropagatedChildren;
    int componentIndex;

    MemberPathType pathType;
    reduced_member_id nextPathMember;
    bool nextPathMemberIsParent;
    spqr_arc pathSourceArc;
    spqr_arc pathTargetArc;
} MATRECColReducedMember;

typedef struct {
    int rootDepth;
    reduced_member_id root;

    reduced_member_id pathEndMembers[2];
    int numPathEndMembers;
} MATRECColReducedComponent;

typedef struct {
    reduced_member_id reducedMember;
    reduced_member_id rootDepthMinimizer;
} MemberInfo;

typedef struct {
    spqr_member member;
} CreateReducedMembersCallstack;

struct MATRECNetworkColumnAdditionImpl {
    bool remainsNetwork;

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

    PathArcListNode *pathArcs;
    int memPathArcs;
    int numPathArcs;
    path_arc_id firstOverallPathArc;

    int *nodeInPathDegree;
    int *nodeOutPathDegree;
    int memNodePathDegree;

    bool *arcInPath;
    bool *arcInPathReversed;
    int memArcsInPath;

    CreateReducedMembersCallstack * createReducedMembersCallStack;
    int memCreateReducedMembersCallStack;

    MATREC_col newColIndex;

    MATREC_row *newRowArcs;
    bool * newRowArcReversed;
    int memNewRowArcs;
    int numNewRowArcs;

    spqr_arc *decompositionRowArcs;
    bool *decompositionArcReversed;
    int memDecompositionRowArcs;
    int numDecompositionRowArcs;

    spqr_member * leafMembers;
    int numLeafMembers;
    int memLeafMembers;
};

static void cleanupPreviousIteration(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol) {
    assert(dec);
    assert(newCol);

    path_arc_id pathArc = newCol->firstOverallPathArc;
    while (pathArcIsValid(pathArc)) {
        spqr_node head = newCol->pathArcs[pathArc].arcHead;
        spqr_node tail = newCol->pathArcs[pathArc].arcTail;
        if(SPQRnodeIsValid(head)){
            newCol->nodeInPathDegree[head] = 0;
        }
        if(SPQRnodeIsValid(tail)){
            newCol->nodeOutPathDegree[tail] = 0;
        }

        spqr_arc arc = newCol->pathArcs[pathArc].arc;
        if(arc < newCol->memArcsInPath){
            newCol->arcInPath[arc] = false;
            newCol->arcInPathReversed[arc] = false;
        }
        pathArc = newCol->pathArcs[pathArc].nextOverall;
    }
#ifndef NDEBUG
    for (int i = 0; i < newCol->memArcsInPath; ++i) {
        assert(newCol->arcInPath[i] == false);
        assert(newCol->arcInPathReversed[i] == false);
    }

    for (int i = 0; i < newCol->memNodePathDegree; ++i) {
        assert(newCol->nodeInPathDegree[i] == 0);
        assert(newCol->nodeOutPathDegree[i] == 0);
    }
#endif

    newCol->firstOverallPathArc = INVALID_PATH_ARC;
    newCol->numPathArcs = 0;
}

MATREC_ERROR MATRECcreateNetworkColumnAddition(MATREC *env, MATRECNetworkColumnAddition **pNewCol) {
    assert(env);

    MATREC_CALL(MATRECallocBlock(env, pNewCol));
    MATRECNetworkColumnAddition *newCol = *pNewCol;

    newCol->remainsNetwork = false;
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

    newCol->pathArcs = NULL;
    newCol->memPathArcs = 0;
    newCol->numPathArcs = 0;
    newCol->firstOverallPathArc = INVALID_PATH_ARC;

    newCol->nodeInPathDegree = NULL;
    newCol->nodeOutPathDegree = NULL;
    newCol->memNodePathDegree = 0;

    newCol->arcInPath = NULL;
    newCol->arcInPathReversed = NULL;
    newCol->memArcsInPath = 0;

    newCol->createReducedMembersCallStack = NULL;
    newCol->memCreateReducedMembersCallStack = 0;

    newCol->newColIndex = MATREC_INVALID_COL;

    newCol->newRowArcs = NULL;
    newCol->newRowArcReversed = NULL;
    newCol->memNewRowArcs = 0;
    newCol->numNewRowArcs = 0;

    newCol->decompositionRowArcs = NULL;
    newCol->decompositionArcReversed = NULL;
    newCol->memDecompositionRowArcs = 0;
    newCol->numDecompositionRowArcs = 0;

    newCol->leafMembers = NULL;
    newCol->numLeafMembers = 0;
    newCol->memLeafMembers = 0;

    return MATREC_OKAY;
}

void MATRECfreeNetworkColumnAddition(MATREC *env, MATRECNetworkColumnAddition **pNewCol) {
    assert(env);
    MATRECNetworkColumnAddition *newCol = *pNewCol;
    MATRECfreeBlockArray(env, &newCol->decompositionRowArcs);
    MATRECfreeBlockArray(env, &newCol->decompositionArcReversed);
    MATRECfreeBlockArray(env, &newCol->newRowArcs);
    MATRECfreeBlockArray(env, &newCol->newRowArcReversed);
    MATRECfreeBlockArray(env, &newCol->createReducedMembersCallStack);
    MATRECfreeBlockArray(env, &newCol->arcInPath);
    MATRECfreeBlockArray(env, &newCol->arcInPathReversed);
    MATRECfreeBlockArray(env, &newCol->nodeInPathDegree);
    MATRECfreeBlockArray(env, &newCol->nodeOutPathDegree);
    MATRECfreeBlockArray(env, &newCol->pathArcs);
    MATRECfreeBlockArray(env, &newCol->childrenStorage);
    MATRECfreeBlockArray(env, &newCol->memberInformation);
    MATRECfreeBlockArray(env, &newCol->reducedComponents);
    MATRECfreeBlockArray(env, &newCol->reducedMembers);

    MATRECfreeBlock(env, pNewCol);
}


static reduced_member_id createReducedMembersToRoot(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol, const spqr_member firstMember ){
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

            reducedMemberData->type = REDUCEDMEMBER_TYPE_UNASSIGNED;
            reducedMemberData->numPropagatedChildren = 0;
            reducedMemberData->firstPathArc = INVALID_PATH_ARC;
            reducedMemberData->numPathArcs = 0;
            reducedMemberData->rigidPathStart = SPQR_INVALID_NODE;
            reducedMemberData->rigidPathEnd = SPQR_INVALID_NODE;

            reducedMemberData->componentIndex = -1;
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
                reducedMemberData->componentIndex = newCol->numReducedComponents;

                assert(newCol->numReducedComponents < newCol->memReducedComponents);
                newCol->reducedComponents[newCol->numReducedComponents].root = reducedMember;
                newCol->reducedComponents[newCol->numReducedComponents].numPathEndMembers = 0;
                newCol->reducedComponents[newCol->numReducedComponents].pathEndMembers[0] = INVALID_REDUCED_MEMBER;
                newCol->reducedComponents[newCol->numReducedComponents].pathEndMembers[1] = INVALID_REDUCED_MEMBER;
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
            //ensure that all newly created reduced members are pointing to the correct component
            assert(parentReducedMemberData->componentIndex >= 0);
            reducedMemberData->componentIndex = parentReducedMemberData->componentIndex;

            newCol->reducedMembers[parentReducedMember].numChildren++;
        }

    }

    reduced_member_id returnedMember = newCol->memberInformation[callstack[0].member].reducedMember;
    return returnedMember;
}

static MATREC_ERROR constructReducedDecomposition(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol) {
    assert(dec);
    assert(newCol);
#ifndef NDEBUG
    for (int i = 0; i < newCol->memMemberInformation; ++i) {
        assert(reducedMemberIsInvalid(newCol->memberInformation[i].reducedMember));
    }
#endif
    newCol->numReducedComponents = 0;
    newCol->numReducedMembers = 0;
    if (newCol->numDecompositionRowArcs == 0) { //Early return in case the reduced decomposition will be empty
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
    for (int i = 0; i < newCol->numDecompositionRowArcs; ++i) {
        assert(i < newCol->memDecompositionRowArcs);
        spqr_arc arc = newCol->decompositionRowArcs[i];
        spqr_member arcMember = findArcMember(dec, arc);
        reduced_member_id reducedMember = createReducedMembersToRoot(dec, newCol, arcMember);
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

static void cleanUpMemberInformation(MATRECNetworkColumnAddition * newCol){
    //This loop is at the end as memberInformation is also used to assign the cut arcs during propagation
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

static void createPathArc(
        MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition * newCol,
        const spqr_arc arc, const reduced_member_id reducedMember, bool reversed){
    assert(dec);
    assert(newCol);

    path_arc_id path_arc = newCol->numPathArcs;
    PathArcListNode * listNode = &newCol->pathArcs[path_arc];
    listNode->arc = arc;

    listNode->nextMember = newCol->reducedMembers[reducedMember].firstPathArc;
    newCol->reducedMembers[reducedMember].firstPathArc = path_arc;
    newCol->reducedMembers[reducedMember].numPathArcs += 1;

    listNode->nextOverall = newCol->firstOverallPathArc;
    newCol->firstOverallPathArc = path_arc;

    ++newCol->numPathArcs;
    assert(newCol->numPathArcs <= newCol->memPathArcs);

    assert(arc < newCol->memArcsInPath);
    newCol->arcInPath[arc] = true;
    newCol->arcInPathReversed[arc] = reversed;
    assert(memberIsRepresentative(dec,newCol->reducedMembers[reducedMember].member));
    if(getMemberType(dec,newCol->reducedMembers[reducedMember].member) == SPQR_MEMBERTYPE_RIGID){

        listNode->arcHead = findEffectiveArcHead(dec,arc);
        listNode->arcTail = findEffectiveArcTail(dec,arc);
        if(reversed){
            swap_ints(&listNode->arcHead,&listNode->arcTail);
        }
        assert(SPQRnodeIsValid(listNode->arcHead) && SPQRnodeIsValid(listNode->arcTail));
        assert(listNode->arcHead < newCol->memNodePathDegree && listNode->arcTail < newCol->memNodePathDegree);
        ++newCol->nodeInPathDegree[listNode->arcHead];
        ++newCol->nodeOutPathDegree[listNode->arcTail];
    }else{
        listNode->arcHead = SPQR_INVALID_NODE;
        listNode->arcTail = SPQR_INVALID_NODE;
    }
    listNode->reversed = reversed;

}

static MATREC_ERROR createPathArcs(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol){
    int maxNumPathArcs = newCol->numDecompositionRowArcs + getNumMembers(dec);
    if(newCol->memPathArcs < maxNumPathArcs){
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->pathArcs,(size_t) maxNumPathArcs)); //TODO: fix reallocation strategy
        newCol->memPathArcs = maxNumPathArcs;
    }
    int maxPathArcIndex = largestArcID(dec);
    if(newCol->memArcsInPath < maxPathArcIndex){
        int newSize = maxPathArcIndex;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->arcInPath,(size_t) newSize));//TODO: fix reallocation strategy
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->arcInPathReversed,(size_t) newSize));//TODO: fix reallocation strategy

        for (int i = newCol->memArcsInPath; i < newSize; ++i) {
            newCol->arcInPath[i] = false;
            newCol->arcInPathReversed[i] = false;
        }
        newCol->memArcsInPath = newSize;
    }
    int maxNumNodes = largestNodeID(dec);
    if(newCol->memNodePathDegree < maxNumNodes){
        int newSize = maxNumNodes;
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->nodeInPathDegree,(size_t) newSize));
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newCol->nodeOutPathDegree,(size_t) newSize));
        for (int i = newCol->memNodePathDegree; i < newSize; ++i) {
            newCol->nodeInPathDegree[i] = 0;
            newCol->nodeOutPathDegree[i] = 0;
        }
        newCol->memNodePathDegree = newSize;
    }
    for (int i = 0; i < newCol->numDecompositionRowArcs; ++i) {
        spqr_arc arc = newCol->decompositionRowArcs[i];
        spqr_member member = findArcMember(dec, arc);
        reduced_member_id reducedMember = newCol->memberInformation[member].reducedMember;
        createPathArc(dec,newCol,arc,reducedMember,newCol->decompositionArcReversed[i]);
    }

    return MATREC_OKAY;
}


/**
 * Saves the information of the current row and partitions it based on whether or not the given columns are
 * already part of the decomposition.
 */
static MATREC_ERROR
newColUpdateColInformation(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol, MATREC_col column,
                           const MATREC_row * nonzeroRows, const double * nonzeroValues, size_t numNonzeros) {
    newCol->newColIndex = column;

    newCol->numDecompositionRowArcs = 0;
    newCol->numNewRowArcs = 0;

    for (size_t i = 0; i < numNonzeros; ++i) {
        spqr_arc rowArc = getDecompositionRowArc(dec, nonzeroRows[i]);
        bool reversed = nonzeroValues[i] < 0.0;
        if (SPQRarcIsValid(rowArc)) { //If the arc is the current decomposition: save it in the array
            if (newCol->numDecompositionRowArcs == newCol->memDecompositionRowArcs) {
                int newNumArcs = newCol->memDecompositionRowArcs == 0 ? 8 : 2 *
                                                                              newCol->memDecompositionRowArcs; //TODO: make reallocation numbers more consistent with rest?
                newCol->memDecompositionRowArcs = newNumArcs;
                MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->decompositionRowArcs,
                                                (size_t) newCol->memDecompositionRowArcs));
                MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->decompositionArcReversed,
                                                (size_t) newCol->memDecompositionRowArcs));
            }
            newCol->decompositionRowArcs[newCol->numDecompositionRowArcs] = rowArc;
            newCol->decompositionArcReversed[newCol->numDecompositionRowArcs] = reversed;
            ++newCol->numDecompositionRowArcs;
        } else {
            //Not in the decomposition: add it to the set of arcs which are newly added with this row.
            if (newCol->numNewRowArcs == newCol->memNewRowArcs) {
                int newNumArcs = newCol->memNewRowArcs == 0 ? 8 : 2 *
                                                                    newCol->memNewRowArcs; //TODO: make reallocation numbers more consistent with rest?
                newCol->memNewRowArcs = newNumArcs;
                MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->newRowArcs,
                                                (size_t) newCol->memNewRowArcs));
                MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->newRowArcReversed,
                                                (size_t) newCol->memNewRowArcs));
            }
            newCol->newRowArcs[newCol->numNewRowArcs] = nonzeroRows[i];
            newCol->newRowArcReversed[newCol->numNewRowArcs] = reversed;
            newCol->numNewRowArcs++;
        }
    }

    return MATREC_OKAY;
}

static MATREC_ERROR computeLeafMembers(const MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol) {
    if (newCol->numReducedMembers > newCol->memLeafMembers) {
        newCol->memLeafMembers = max(newCol->numReducedMembers, 2 * newCol->memLeafMembers);
        MATREC_CALL(MATRECreallocBlockArray(dec->env, &newCol->leafMembers, (size_t) newCol->memLeafMembers));
    }
    newCol->numLeafMembers = 0;

    for (reduced_member_id reducedMember = 0; reducedMember < newCol->numReducedMembers; ++reducedMember) {
        if (newCol->reducedMembers[reducedMember].numChildren == 0) {
            newCol->leafMembers[newCol->numLeafMembers] = reducedMember;
            ++newCol->numLeafMembers;
        }
    }
    return MATREC_OKAY;
}

static void determineRigidPath(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol, MATRECColReducedMember * redMem){
    assert(dec);
    assert(newCol);
    assert(redMem);

    bool isValidPath = true;
    redMem->rigidPathStart = SPQR_INVALID_NODE;
    redMem->rigidPathEnd = SPQR_INVALID_NODE;
    for(path_arc_id pathArc = redMem->firstPathArc; pathArcIsValid(pathArc); pathArc = newCol->pathArcs[pathArc].nextMember){
        spqr_node head = newCol->pathArcs[pathArc].arcHead;
        spqr_node tail = newCol->pathArcs[pathArc].arcTail;
        assert(nodeIsRepresentative(dec,head) && nodeIsRepresentative(dec,tail));

        if(newCol->nodeInPathDegree[head] > 1 || newCol->nodeOutPathDegree[tail] > 1){
            //not network -> stop
            isValidPath = false;
            break;
        }
        if(newCol->nodeOutPathDegree[head] == 0){
            //found end node
            //If this is the second, stop
            if(SPQRnodeIsValid(redMem->rigidPathEnd)){
                isValidPath = false;
                break;
            }
            redMem->rigidPathEnd = head;
        }
        if(newCol->nodeInPathDegree[tail] == 0){
            //Found start node.
            //If this is the second, stop.
            if(SPQRnodeIsValid(redMem->rigidPathStart)){
                isValidPath = false;
                break;
            }
            redMem->rigidPathStart = tail;
        }
    }
    //assert that both a start and end node have been found
    assert(!isValidPath || (SPQRnodeIsValid(redMem->rigidPathStart) && SPQRnodeIsValid(redMem->rigidPathEnd)));
    if(!isValidPath){
        redMem->rigidPathStart = SPQR_INVALID_NODE;
        redMem->rigidPathEnd = SPQR_INVALID_NODE;
        redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
        newCol->remainsNetwork = false;
    }

}
static void determineSingleRigidType(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol,
                                     reduced_member_id reducedMember){
    assert(dec);
    assert(newCol);
    assert(reducedMemberIsValid(reducedMember));
    MATRECColReducedMember * redMem = &newCol->reducedMembers[reducedMember];
    assert(pathArcIsValid(redMem->firstPathArc));
    determineRigidPath(dec,newCol,redMem);
    if(redMem->type != REDUCEDMEMBER_TYPE_NOT_NETWORK){
        redMem->type = REDUCEDMEMBER_TYPE_MERGED;
    }
}
//TODO: type seems somewhat duplicate
static void determineSingleComponentType(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol, reduced_member_id reducedMember){
    assert(dec);
    assert(newCol);

    int numNonPropagatedAdjacent = newCol->reducedMembers[reducedMember].numChildren-newCol->reducedMembers[reducedMember].numPropagatedChildren;
    if(reducedMemberIsValid(newCol->reducedMembers[reducedMember].parent) &&
            newCol->reducedMembers[newCol->reducedMembers[reducedMember].parent].type != REDUCEDMEMBER_TYPE_CYCLE){
        ++numNonPropagatedAdjacent;
    }

    if(numNonPropagatedAdjacent > 2){
        newCol->reducedMembers[reducedMember].type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
        newCol->remainsNetwork = false;
        return;
    }


    spqr_member member = findMember(dec, newCol->reducedMembers[reducedMember].member); //TODO: find necessary?
    SPQRMemberType type = getMemberType(dec, member);
    switch(type){
        case SPQR_MEMBERTYPE_RIGID:
        {
            determineSingleRigidType(dec,newCol,reducedMember);
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:
        {
            MATRECColReducedMember *redMem =&newCol->reducedMembers[reducedMember];
            assert(pathArcIsValid(redMem->firstPathArc));
            bool pathForwards = newCol->pathArcs[redMem->firstPathArc].reversed ==
                    arcIsReversedNonRigid(dec,newCol->pathArcs[redMem->firstPathArc].arc);
            redMem->pathBackwards = !pathForwards;
            redMem->type = REDUCEDMEMBER_TYPE_CYCLE;
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:
        case SPQR_MEMBERTYPE_LOOP:
        {
            MATRECColReducedMember *redMem =&newCol->reducedMembers[reducedMember];
            int countedPathArcs = 0;
            bool good = true;
            bool passesForwards = true;
            for(path_arc_id pathArc = redMem->firstPathArc; pathArcIsValid(pathArc);
                pathArc = newCol->pathArcs[pathArc].nextMember){
                if(countedPathArcs == 0){
                    passesForwards = newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec,newCol->pathArcs[pathArc].arc);
                }else if((newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec,newCol->pathArcs[pathArc].arc)) != passesForwards){
                    good = false;
                    break;
                }
                ++countedPathArcs;
            }
            if(!good){
                redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
                newCol->remainsNetwork = false;
                break;
            }

            redMem->pathBackwards = !passesForwards;
            if (countedPathArcs == getNumMemberArcs(dec, findMember(dec,redMem->member)) -1){
                //Type -> Cycle;
                //Propagate arc
                redMem->type = REDUCEDMEMBER_TYPE_CYCLE;
            }else{
                //Type -> single_end
                redMem->type = REDUCEDMEMBER_TYPE_MERGED;
            }
            break;
        }
        case SPQR_MEMBERTYPE_UNASSIGNED:
        {
            assert(false);
            break;
        }
    }
}


static void determinePathSeriesType(MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition *newCol,
                                    reduced_member_id reducedMember, spqr_member member, MemberPathType previousType,
                                    spqr_arc source, spqr_arc target){
    assert(dec);
    assert(newCol);
    assert(reducedMemberIsValid(reducedMember));
    assert(SPQRmemberIsValid(member) &&memberIsRepresentative(dec,member));
    assert(getMemberType(dec,member) == SPQR_MEMBERTYPE_SERIES);

    MATRECColReducedMember *redMem =&newCol->reducedMembers[reducedMember];
    int countedPathArcs = 0;

    bool good = true;
    bool passesForwards = true;
    for(path_arc_id pathArc = redMem->firstPathArc; pathArcIsValid(pathArc);
        pathArc = newCol->pathArcs[pathArc].nextMember){
        if(countedPathArcs == 0){
            passesForwards = newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec,newCol->pathArcs[pathArc].arc);
        }else if((newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec,newCol->pathArcs[pathArc].arc)) != passesForwards){
            good = false;
            break;
        }
        ++countedPathArcs;
    }
    //If the internal directions of the arcs do not agree, we have no way to realize
    if(!good){
        redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
        newCol->remainsNetwork = false;
        return;
    }
    //If we are in the first skeleton processed, ignore the previous member type
    if(!SPQRarcIsValid(source)){
        assert(countedPathArcs > 0);
        assert(SPQRarcIsValid(target));
        bool firstReversed = arcIsReversedNonRigid(dec,newCol->pathArcs[redMem->firstPathArc].arc);
        bool targetReversed = arcIsReversedNonRigid(dec,target);
        bool reversePath = newCol->pathArcs[redMem->firstPathArc].reversed;

        if((firstReversed == targetReversed ) == reversePath){
            redMem->pathType = INTO_HEAD;
        }else{
            redMem->pathType = OUT_HEAD;
        }
        redMem->pathBackwards = !passesForwards;
        return;
    }
    bool sourceReversed = arcIsReversedNonRigid(dec,source);
    if(countedPathArcs > 0){
        bool isIntoHeadOrOutTail = isInto(previousType) == isHead(previousType);
        bool isGood = isIntoHeadOrOutTail == (sourceReversed == passesForwards); //TODO: check if this is correct
        if(!isGood){
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
            return;
        }
    }
    redMem->pathBackwards = !passesForwards;
    if(SPQRarcIsValid(target)){
        bool targetReversed = arcIsReversedNonRigid(dec,target);
        bool inSameDirection = sourceReversed == targetReversed;

        MemberPathType currentType;
        switch(previousType){
            case INTO_HEAD:{
                currentType = inSameDirection ? INTO_TAIL : INTO_HEAD;
                break;
            }
            case INTO_TAIL:
            {
                currentType = inSameDirection ? INTO_HEAD : INTO_TAIL;
                break;
            }
            case OUT_HEAD:{
                currentType = inSameDirection ? OUT_TAIL : OUT_HEAD;
                break;
            }
            case OUT_TAIL:
            {
                currentType = inSameDirection ? OUT_HEAD : OUT_TAIL;
                break;
            }
            default:{
                assert(false);
            }
        }
        redMem->pathType = currentType;
        return;
    }
    //If we are in the last skeleton, we only have a source, so nothing further to do
    assert(countedPathArcs > 0);
    //Strictly speaking below are no-ops within the algorithm, but first we focus on not getting any bugs (TODO)
    if(isInto(previousType)){
        redMem->pathType = INTO_HEAD;
    }else{
        redMem->pathType = OUT_HEAD;
    }
}
static void determinePathParallelType(MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition *newCol,
                                      reduced_member_id reducedMember, spqr_member member, MemberPathType previousType,
                                      spqr_arc source, spqr_arc target) {
    assert(dec);
    assert(newCol);
    assert(reducedMemberIsValid(reducedMember));
    assert(SPQRmemberIsValid(member) && memberIsRepresentative(dec, member));
    assert(getMemberType(dec, member) == SPQR_MEMBERTYPE_PARALLEL);

    //Parallel members must always be of degree two; if they are a leaf, they must contain an edge, which could have been propagated
    assert(SPQRarcIsValid(source) && SPQRarcIsValid(target));
    bool sourceReversed = arcIsReversedNonRigid(dec,source);
    bool targetReversed = arcIsReversedNonRigid(dec,target);
    bool inSameDirection = sourceReversed == targetReversed;

    MATRECColReducedMember *redMem =&newCol->reducedMembers[reducedMember];

    path_arc_id pathArc = redMem->firstPathArc;

    bool arcPresent = false;
    if(pathArcIsValid(pathArc)){
        bool pathArcReversed = newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec,newCol->pathArcs[pathArc].arc);
        bool intoHeadOrOutTail = isInto(previousType) == isHead(previousType);
        bool good = intoHeadOrOutTail == (pathArcReversed != sourceReversed);
        if(!good){
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
            return;
        }
        arcPresent = true;
    }

    bool swapHeadTail = arcPresent != inSameDirection;
    MemberPathType currentType;
    switch(previousType){
        case INTO_HEAD:{
            currentType = swapHeadTail ? INTO_HEAD : INTO_TAIL;
            break;
        }
        case INTO_TAIL:{
            currentType = swapHeadTail ? INTO_TAIL : INTO_HEAD;
            break;
        }
        case OUT_HEAD:{
            currentType = swapHeadTail ? OUT_HEAD  : OUT_TAIL;
            break;
        }
        case OUT_TAIL:{
            currentType = swapHeadTail ? OUT_TAIL : OUT_HEAD;
            break;
        }
    }
    redMem->pathType = currentType;
    //TODO: set path reversed?
}
static void determinePathRigidType(MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition *newCol,
                                   reduced_member_id reducedMember, MemberPathType previousType,
                                   spqr_arc source, spqr_arc target) {

    MATRECColReducedMember *redMem = &newCol->reducedMembers[reducedMember];
    if(pathArcIsInvalid(redMem->firstPathArc)){
        assert(SPQRarcIsValid(source));
        assert(SPQRarcIsValid(target));
        //In this case, we need to check if the source and target are adjacent in any node

        spqr_node sourceTail = findEffectiveArcTail(dec, source);
        spqr_node sourceHead = findEffectiveArcHead(dec, source);
        spqr_node targetTail = findEffectiveArcTail(dec, target);
        spqr_node targetHead = findEffectiveArcHead(dec, target);
        bool sourceHeadIsTargetHead = sourceHead == targetHead;
        bool sourceTailIsTargetHead = sourceTail == targetHead;
        bool sourceHeadIsTargetTail = sourceHead == targetTail;
        bool sourceTailIsTargetTail = sourceTail == targetTail;

        if(!(sourceHeadIsTargetHead || sourceTailIsTargetHead || sourceHeadIsTargetTail || sourceTailIsTargetTail)){
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
            return;
        }
        assert((sourceHeadIsTargetHead ? 1 : 0) + (sourceHeadIsTargetTail ? 1 : 0)
        + (sourceTailIsTargetHead ? 1 : 0) +  (sourceTailIsTargetTail ? 1 : 0) == 1 );
        //assert simplicity; they can be adjacent in only exactly one node
        bool isSourceHead = sourceHeadIsTargetHead || sourceHeadIsTargetTail;
        bool isTargetTail = sourceHeadIsTargetTail || sourceTailIsTargetTail;
        if (isHead(previousType) == isSourceHead){
            //no need to reflect
            redMem->reverseArcs = false;
            if(isInto(previousType)){
                if(isTargetTail){
                    redMem->pathType = INTO_TAIL;
                }else{
                    redMem->pathType = INTO_HEAD;
                }
            }else{
                if(isTargetTail){
                    redMem->pathType = OUT_TAIL;
                }else{
                    redMem->pathType = OUT_HEAD;
                }
            }
        }else{
            redMem->reverseArcs = true;
            //because of the reversal, all the heads/tails are switched below
            if(isInto(previousType)){
                if(isTargetTail){
                    redMem->pathType = INTO_HEAD;
                }else{
                    redMem->pathType = INTO_TAIL;
                }
            }else{
                if(isTargetTail){
                    redMem->pathType = OUT_HEAD;
                }else{
                    redMem->pathType = OUT_TAIL;
                }
            }
        }
        assert(isInto(redMem->pathType) == isInto(previousType));
        return;
    }
    determineRigidPath(dec, newCol, redMem);
    if (redMem->type == REDUCEDMEMBER_TYPE_NOT_NETWORK) {
        return;
    }
    if(SPQRarcIsInvalid(source)){
        assert(SPQRarcIsValid(target));
        spqr_node targetTail = findEffectiveArcTail(dec, target);
        spqr_node targetHead = findEffectiveArcHead(dec, target);
        redMem->reverseArcs = false;
        if(redMem->rigidPathEnd == targetHead){
            redMem->pathType = INTO_HEAD;
        }else if(redMem->rigidPathEnd == targetTail){
            redMem->pathType = INTO_TAIL;
        }else if(redMem->rigidPathStart == targetHead){
            redMem->pathType = OUT_HEAD;
        }else if(redMem->rigidPathStart == targetTail){
            redMem->pathType = OUT_TAIL;
        }else {
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
        }
        return;
    }
    assert(SPQRarcIsValid(source));
    spqr_node sourceTail = findEffectiveArcTail(dec, source);
    spqr_node sourceHead = findEffectiveArcHead(dec, source);

    bool startsAtHead = sourceHead == redMem->rigidPathStart;
    bool endsAtTail = sourceTail == redMem->rigidPathEnd;
    bool startsAtTail = sourceTail == redMem->rigidPathStart;
    bool endsAtHead = sourceHead == redMem->rigidPathEnd;

    bool isIntoHeadOrOutTail = isInto(previousType) == isHead(previousType);
    if(isIntoHeadOrOutTail){ //into head or outTail
        //Check if path starts at head or ends at tail
        if(!startsAtHead && !endsAtTail){
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
            return;
        }
        assert(startsAtHead || endsAtTail); //both can hold; they can form cycle but other components can not be reduced
//        redMem->reverseArcs = isInto(previousType) != startsAtHead; //Reverse only if there is no path starting at head
    }else{ //Into tail or outHead
        //Check if path starts at tail or ends at head
        if(!startsAtTail && !endsAtHead){
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
            return;
        }
        assert(startsAtTail || endsAtHead); // both can hold; they can form cycle but other components can not be reduced
//        redMem->reverseArcs = isInto(previousType) != startsAtTail; //Reverse only if there is no path starting at tail
    }

    if(SPQRarcIsValid(target)){
        spqr_node targetTail = findEffectiveArcTail(dec, target);
        spqr_node targetHead = findEffectiveArcHead(dec, target);

        //Check if they are not parallel; (below logic relies on this fact)
        assert(!((targetHead == sourceHead && targetTail == sourceTail) || (targetHead == sourceTail && targetTail == sourceHead)));

        bool startsAtTargetHead = redMem->rigidPathStart == targetHead;
        bool startsAtTargetTail = redMem->rigidPathStart == targetTail;
        bool endsAtTargetHead = redMem->rigidPathEnd == targetHead;
        bool endsAtTargetTail = redMem->rigidPathEnd == targetTail;

        if(!(startsAtTargetHead || startsAtTargetTail || endsAtTargetHead || endsAtTargetTail)){
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
            return;
        }
        bool outReverse = false;
        bool outHead = false;

        if(isInto(previousType) == isHead(previousType)) {
            if (startsAtHead && endsAtTail) {
                outReverse = (startsAtTargetHead || startsAtTargetTail) == isInto(previousType);
            } else if(startsAtHead) {
                outReverse = !isInto(previousType);
            } else {
                assert(endsAtTail);
                outReverse = isInto(previousType);
            }
        }else{
            if (startsAtTail && endsAtHead) {
                outReverse = (startsAtTargetHead || startsAtTargetTail) == isInto(previousType);
            } else if (startsAtTail) {
                outReverse = !isInto(previousType);
            } else {
                assert(endsAtHead);
                outReverse = isInto(previousType);
            }
        }

        //TODO: these can probably be simplified significantly, but this might pose risk of introducing incorrect assumptions
        bool isBad = false;
        if(isInto(previousType) == isHead(previousType)){
            if(startsAtHead && endsAtTail){
                outHead = (startsAtTargetTail || endsAtTargetHead) == isInto(previousType);
            }else if(startsAtHead){
                if(endsAtTargetHead){
                    outHead = isInto(previousType);
                }else if (endsAtTargetTail){
                    outHead = !isInto(previousType);
                }else{
                    isBad = true;
                }
            }else{
                assert(endsAtTail);
                if(startsAtTargetTail){
                    outHead = isInto(previousType);
                }else if (startsAtTargetHead){
                    outHead = !isInto(previousType);
                }else{
                    isBad = true;
                }
            }
        }else{
            if(startsAtTail && endsAtHead){
                outHead = (startsAtTargetTail || endsAtTargetHead) == isInto(previousType);
            }else if(startsAtTail){
                if(endsAtTargetHead){
                    outHead = isInto(previousType);
                }else if(endsAtTargetTail){
                    outHead = !isInto(previousType);
                }else{
                    isBad = true;
                }
            }else{
                assert(endsAtHead);
                if(startsAtTargetTail){
                    outHead = isInto(previousType);
                }else if (startsAtTargetHead){
                    outHead = !isInto(previousType);
                }else{
                    isBad = true;
                }
            }
        }
        if(isBad){
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = false;
            return;
        }

        redMem->reverseArcs = outReverse;
        if(isInto(previousType)){
            redMem->pathType = outHead ? INTO_HEAD : INTO_TAIL;
        }else{
            redMem->pathType = outHead ? OUT_HEAD : OUT_TAIL;
        }
        return;
    }

    //TODO: is this simplifyable?
    if(isInto(previousType) == isHead(previousType)){
        redMem->reverseArcs = startsAtHead != isInto(previousType);
    }else{
        redMem->reverseArcs = startsAtTail != isInto(previousType);
    }
    //last member of the path. Since we already checked the source,
    //Below is technically no op (TODO)
    if(isInto(previousType)){
        redMem->pathType = INTO_HEAD;
    }else{
        redMem->pathType = OUT_HEAD;
    }
}
static void determinePathMemberType(MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition *newCol,
                                    reduced_member_id reducedMember, spqr_member member, MemberPathType previousType,
                                    spqr_arc source, spqr_arc target){
    newCol->reducedMembers[reducedMember].pathSourceArc = source;
    newCol->reducedMembers[reducedMember].pathTargetArc = target;
    //Check if the marked edges with the given signs
    //form a (reverse) directed path from one of the source's end nodes to one of the target's end nodes
    switch(getMemberType(dec,member)){
        case SPQR_MEMBERTYPE_RIGID:{
            determinePathRigidType(dec,newCol,reducedMember,previousType,source,target);
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:{
            determinePathParallelType(dec,newCol,reducedMember,member,previousType,source,target);
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:{
            determinePathSeriesType(dec,newCol,reducedMember,member,previousType,source,target);
            break;
        }
        case SPQR_MEMBERTYPE_LOOP:
        case SPQR_MEMBERTYPE_UNASSIGNED:
        {
            assert(false);
            //In release,

            newCol->remainsNetwork = false;
            break;
        }
    }

}

static void determinePathTypes(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol, MATRECColReducedComponent * component){
    assert(dec);
    assert(newCol);
    assert(component);
    assert(component->numPathEndMembers == 2);

    assert(component->pathEndMembers[0] != component->root);

    //We check the path by going from end to end. We start at the leaf and process every component,
    //walking down until we hit the root.
    //Then, we walk up from the root and process each component in the same manner.
    reduced_member_id reducedStart = component->pathEndMembers[0];
    spqr_arc toPrevious = SPQR_INVALID_ARC;
    spqr_arc toNext = SPQR_INVALID_ARC;
    MemberPathType previousType = INTO_HEAD; //Arbitrary, is ignored in the first call

    spqr_member member = newCol->reducedMembers[reducedStart].member;
    reduced_member_id reducedMember = reducedStart;
    reduced_member_id previousReducedMember = reducedStart;
    while(reducedMember != component->root){
        toNext = markerToParent(dec,member);
        determinePathMemberType(dec,newCol,reducedMember,member,previousType,toPrevious,toNext);
        if(!newCol->remainsNetwork){
            return;
        }
        previousType = newCol->reducedMembers[reducedMember].pathType;
        toPrevious = markerOfParent(dec,member);
        member = findMemberParent(dec,newCol->reducedMembers[reducedMember].member);
        previousReducedMember = reducedMember;
        reducedMember = newCol->memberInformation[member].reducedMember;
        newCol->reducedMembers[previousReducedMember].nextPathMember = reducedMember;
        newCol->reducedMembers[previousReducedMember].nextPathMemberIsParent = true;
    }

    while(reducedMember != component->pathEndMembers[1]){
        //Search the (other) child node
        reduced_member_id child = INVALID_REDUCED_MEMBER;
        //a bit ugly linear search, but not a problem for time complexity
        for (children_idx i = newCol->reducedMembers[reducedMember].firstChild;
        i < newCol->reducedMembers[reducedMember].firstChild + newCol->reducedMembers[reducedMember].numChildren; ++i) {
            reduced_member_id childReduced = newCol->childrenStorage[i];
            if(newCol->reducedMembers[childReduced].type != REDUCEDMEMBER_TYPE_CYCLE &&
               childReduced != previousReducedMember){
                child = childReduced;
                break;
            }
        }
        assert(reducedMemberIsValid(child));

        spqr_member childMember = newCol->reducedMembers[child].member;
        toNext = markerOfParent(dec,childMember);

        determinePathMemberType(dec,newCol,reducedMember,member,previousType,toPrevious,toNext);
        if(!newCol->remainsNetwork){
            return;
        }
        previousType = newCol->reducedMembers[reducedMember].pathType;
        toPrevious = markerToParent(dec,childMember);
        member = childMember;
        previousReducedMember = reducedMember;
        reducedMember = child;
        newCol->reducedMembers[previousReducedMember].nextPathMember = reducedMember;
        newCol->reducedMembers[previousReducedMember].nextPathMemberIsParent = false;
    }
    //The last iteration is not performed by the loops above.
    //We explicitly set the target arc to invalid in order to indicate that this is the last iteration.
    toNext = SPQR_INVALID_ARC;
    determinePathMemberType(dec,newCol,reducedMember,member,previousType,toPrevious,toNext);
    newCol->reducedMembers[reducedMember].nextPathMember = INVALID_REDUCED_MEMBER;
    //since we return anyways, no need to check newCol->remainsNetwork explicitly
}
static void checkRigidLeaf(MATRECNetworkDecomposition * dec,MATRECNetworkColumnAddition * newCol, reduced_member_id leaf,
                           spqr_arc toParent, reduced_member_id parent, spqr_arc toChild){
    MATRECColReducedMember * leafMember = &newCol->reducedMembers[leaf];
    determineRigidPath(dec,newCol,leafMember);
    if(leafMember->type == REDUCEDMEMBER_TYPE_NOT_NETWORK){
        return;
    }
    spqr_node targetHead = findEffectiveArcHead(dec,toParent);
    spqr_node targetTail = findEffectiveArcTail(dec,toParent);
    bool matches = leafMember->rigidPathStart == targetTail && leafMember->rigidPathEnd == targetHead;
    bool opposite = leafMember->rigidPathStart == targetHead && leafMember->rigidPathEnd == targetTail;
    if(matches || opposite) {
        leafMember->type = REDUCEDMEMBER_TYPE_CYCLE;
        createPathArc(dec,newCol,toChild,parent, opposite);
        return;
    }
    leafMember->type = REDUCEDMEMBER_TYPE_MERGED;
}
static ReducedMemberType checkLeaf(MATRECNetworkDecomposition * dec,MATRECNetworkColumnAddition * newCol, reduced_member_id leaf,
                      spqr_arc toParent, reduced_member_id parent, spqr_arc toChild){
    assert(dec);
    assert(newCol);
    assert(SPQRarcIsValid(toParent));
    assert(SPQRarcIsValid(toChild));
    assert(!arcIsTree(dec,toParent));
    assert(reducedMemberIsValid(leaf));
    assert(reducedMemberIsValid(parent));

    switch(getMemberType(dec,newCol->reducedMembers[leaf].member)){
        case SPQR_MEMBERTYPE_RIGID:
        {
            checkRigidLeaf(dec,newCol,leaf,toParent,parent,toChild);
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:
        {
            MATRECColReducedMember *reducedMember =&newCol->reducedMembers[leaf];
            assert(pathArcIsValid(reducedMember->firstPathArc));
            reducedMember->type = REDUCEDMEMBER_TYPE_CYCLE;

            bool pathArcReversed = newCol->pathArcs[reducedMember->firstPathArc].reversed;
            bool arcPathArcIsReverse = arcIsReversedNonRigid(dec, newCol->pathArcs[reducedMember->firstPathArc].arc);
            bool parentReversed = arcIsReversedNonRigid(dec, toParent);
            createPathArc(dec, newCol, toChild, parent, (arcPathArcIsReverse == parentReversed) == pathArcReversed);
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:
        case SPQR_MEMBERTYPE_LOOP:
        {
            MATRECColReducedMember *reducedMember =&newCol->reducedMembers[leaf];
            int countedPathArcs = 0;
            bool good = true;
            bool passesForwards = true;
            for(path_arc_id pathArc = reducedMember->firstPathArc; pathArcIsValid(pathArc);
                pathArc = newCol->pathArcs[pathArc].nextMember){
                if(countedPathArcs == 0){
                    passesForwards = newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec,newCol->pathArcs[pathArc].arc);
                }else if((newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec,newCol->pathArcs[pathArc].arc)) != passesForwards){
                    good = false;
                    break;
                }
                ++countedPathArcs;
            }
            if(!good){
                reducedMember->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
                newCol->remainsNetwork = false;
                break;
            }

            reducedMember->pathBackwards = !passesForwards;
            if (countedPathArcs == getNumMemberArcs(dec, findMember(dec,reducedMember->member)) -1){
                //Type -> Cycle;
                //Propagate arc
                reducedMember->type = REDUCEDMEMBER_TYPE_CYCLE;

                bool firstArcReversed = arcIsReversedNonRigid(dec,newCol->pathArcs[reducedMember->firstPathArc].arc);
                bool firstArcInPathReverse = newCol->pathArcs[reducedMember->firstPathArc].reversed;
                bool parentReversed = arcIsReversedNonRigid(dec,toParent);
                createPathArc(dec,newCol,toChild,parent,(parentReversed == firstArcReversed) != firstArcInPathReverse);

            }else{
                //Type -> single_end
                reducedMember->type = REDUCEDMEMBER_TYPE_MERGED;
            }

            break;
        }
        case SPQR_MEMBERTYPE_UNASSIGNED:
        {
            assert(false);
            break;
        }
    }
    return newCol->reducedMembers[leaf].type;
}

static void propagateCycles(MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition * newCol){
    assert(dec);
    assert(newCol);
    int leafArrayIndex = 0;

    while(leafArrayIndex != newCol->numLeafMembers){
        reduced_member_id leaf = newCol->leafMembers[leafArrayIndex];
        //next is invalid if the member is not in the reduced decomposition.
        reduced_member_id next = newCol->reducedMembers[leaf].parent;
        spqr_arc parentMarker = markerToParent(dec, newCol->reducedMembers[leaf].member);
        if(reducedMemberIsValid(next) && !arcIsTree(dec,parentMarker)){
            assert(reducedMemberIsValid(next));
            assert(SPQRarcIsValid(parentMarker));
            ReducedMemberType type = checkLeaf(dec,newCol,leaf,parentMarker,next,markerOfParent(dec,newCol->reducedMembers[leaf].member));
            if(type == REDUCEDMEMBER_TYPE_CYCLE){
                ++newCol->reducedMembers[next].numPropagatedChildren;
                if(newCol->reducedMembers[next].numPropagatedChildren == newCol->reducedMembers[next].numChildren){
                    newCol->leafMembers[leafArrayIndex] = next;
                }else{
                    ++leafArrayIndex;
                }
            }else if(type == REDUCEDMEMBER_TYPE_NOT_NETWORK){
                return;
            }else{
                assert(type == REDUCEDMEMBER_TYPE_MERGED);
                ++leafArrayIndex;
                int component = newCol->reducedMembers[leaf].componentIndex;
                if(newCol->reducedComponents[component].numPathEndMembers >= 2){
                    newCol->remainsNetwork = false;
                    return;
                }
                assert(newCol->reducedComponents[component].numPathEndMembers < 2);
                newCol->reducedComponents[component].pathEndMembers[newCol->reducedComponents[component].numPathEndMembers] = leaf;
                ++newCol->reducedComponents[component].numPathEndMembers;
            }
        }else{
            ++leafArrayIndex;
            int component = newCol->reducedMembers[leaf].componentIndex;
            if(newCol->reducedComponents[component].numPathEndMembers >= 2){
                newCol->remainsNetwork = false;
                return;
            }
            assert(newCol->reducedComponents[component].numPathEndMembers < 2);
            newCol->reducedComponents[component].pathEndMembers[newCol->reducedComponents[component].numPathEndMembers] = leaf;
            ++newCol->reducedComponents[component].numPathEndMembers;
        }

    }

    for (int j = 0; j < newCol->numReducedComponents; ++j) {
        //The reduced root might be a leaf as well: we propagate it last
        reduced_member_id root = newCol->reducedComponents[j].root;

        while(true){
            if(newCol->reducedMembers[root].numPropagatedChildren != newCol->reducedMembers[root].numChildren -1) {
                break;
            }
            //TODO: bit ugly, have to do a linear search for the child
            reduced_member_id child = INVALID_REDUCED_MEMBER;
            spqr_arc markerToChild = SPQR_INVALID_ARC;
            for (children_idx i = newCol->reducedMembers[root].firstChild;
                 i < newCol->reducedMembers[root].firstChild + newCol->reducedMembers[root].numChildren; ++i) {
                reduced_member_id childReduced = newCol->childrenStorage[i];
                if (newCol->reducedMembers[childReduced].type != REDUCEDMEMBER_TYPE_CYCLE) {
                    child = childReduced;
                    markerToChild = markerOfParent(dec, newCol->reducedMembers[child].member);
                    break;
                }
            }
            assert(SPQRmemberIsValid(newCol->reducedMembers[child].member));
            assert(SPQRarcIsValid(markerToChild));
            if (!arcIsTree(dec, markerToChild)) {
                ReducedMemberType type = checkLeaf(dec, newCol, root, markerToChild, child,
                                                   markerToParent(dec, newCol->reducedMembers[child].member));
                if (type == REDUCEDMEMBER_TYPE_CYCLE) {
                    root = child;
                    continue;
                } else if (type == REDUCEDMEMBER_TYPE_NOT_NETWORK) {
                    return;
                }
            }
            //If the root has exactly one neighbour and is not contained, it is also considered a path end member
            int component = newCol->reducedMembers[root].componentIndex;
            bool rootPresent = false;
            for (int i = 0; i < newCol->reducedComponents[component].numPathEndMembers; ++i) {
                rootPresent = rootPresent || (newCol->reducedComponents[component].pathEndMembers[i] == root);
            }
            if(!rootPresent){
                if(newCol->reducedComponents[component].numPathEndMembers >= 2){
                    newCol->remainsNetwork = false;
                    return;
                }
                newCol->reducedComponents[component].pathEndMembers[newCol->reducedComponents[component].numPathEndMembers] = root;
                ++newCol->reducedComponents[component].numPathEndMembers;
            }
            break;

        }

        newCol->reducedComponents[j].root = root;
        newCol->reducedMembers[root].parent = INVALID_REDUCED_MEMBER;
    }
}
static void determineComponentTypes(MATRECNetworkDecomposition * dec, MATRECNetworkColumnAddition * newCol, MATRECColReducedComponent * component){
    assert(dec);
    assert(newCol);
    assert(component);

    if(component->numPathEndMembers == 1){
        assert(component->root == component->pathEndMembers[0]);
        determineSingleComponentType(dec,newCol,component->root);
    }else{
        assert(component->numPathEndMembers == 2);
        determinePathTypes(dec,newCol,component);
    }
}

MATREC_ERROR
MATRECNetworkColumnAdditionCheck(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol, MATREC_col column, const MATREC_row * nonzeroRows,
                                 const double * nonzeroValues, size_t numNonzeros) {
    assert(dec);
    assert(newCol);
    assert(numNonzeros == 0 || (nonzeroRows && nonzeroValues));

    newCol->remainsNetwork = true;
    cleanupPreviousIteration(dec, newCol);
    //assert that previous iteration was cleaned up

    //Store call data
    MATREC_CALL(newColUpdateColInformation(dec, newCol, column, nonzeroRows, nonzeroValues, numNonzeros));

    //compute reduced decomposition
    MATREC_CALL(constructReducedDecomposition(dec, newCol));
    //initialize path arcs in reduced decomposition
    MATREC_CALL(createPathArcs(dec,newCol));
    MATREC_CALL(computeLeafMembers(dec,newCol));
    propagateCycles(dec,newCol);
    //determine types
    if(newCol->remainsNetwork){
        for (int i = 0; i < newCol->numReducedComponents; ++i) {
            determineComponentTypes(dec,newCol,&newCol->reducedComponents[i]);
        }
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
    spqr_node head;
    spqr_node tail;
    spqr_arc representative;
    bool reversed;
} NewColInformation;

static NewColInformation emptyNewColInformation(void){
    NewColInformation information;
    information.member = SPQR_INVALID_MEMBER;
    information.head = SPQR_INVALID_NODE;
    information.tail = SPQR_INVALID_NODE;
    information.reversed = false;
    information.representative = SPQR_INVALID_ARC;
    return information;
}
static void setTerminalHead(NewColInformation * info, spqr_node node){
    assert(SPQRnodeIsValid(node));
    assert(SPQRnodeIsInvalid(info->head));
    assert(info);
    info->head = node;
}
static void setTerminalTail(NewColInformation * info, spqr_node node){
    assert(SPQRnodeIsValid(node));
    assert(info);
    assert(SPQRnodeIsInvalid(info->tail));
    info->tail = node;
}
static void setTerminalReversed(NewColInformation *info, bool reversed){
    assert(info);
    info->reversed = reversed;
}
static void setTerminalMember(NewColInformation * info, spqr_member member){
    assert(info);
    info->member = member;
}
static void setTerminalRepresentative(NewColInformation * info, spqr_arc representative){
    assert(info);
    info->representative = representative;
}

static MATREC_ERROR splitParallel(MATRECNetworkDecomposition *dec, spqr_member parallel,
                                spqr_arc arc1, spqr_arc arc2,
                                spqr_member * childParallel){
    assert(dec);
    assert(SPQRarcIsValid(arc1));
    assert(SPQRarcIsValid(arc2));
    assert(SPQRmemberIsValid(parallel));

    bool childContainsTree = arcIsTree(dec,arc1) || arcIsTree(dec,arc2);
    spqr_arc toParent = markerToParent(dec, parallel);
    bool parentMoved = toParent == arc1 || toParent == arc2;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, childParallel));

    moveArcToNewMember(dec,arc1,parallel,*childParallel);
    moveArcToNewMember(dec,arc2,parallel,*childParallel);

    if(parentMoved){
        MATREC_CALL(createMarkerPair(dec,*childParallel,parallel,!childContainsTree,false,false));
    }else{
        MATREC_CALL(createMarkerPair(dec,parallel,*childParallel,childContainsTree,false,false));
    }
    return MATREC_OKAY;
}

static MATREC_ERROR splitSeries(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                              MATRECColReducedMember * reducedMember,
                              spqr_member member, spqr_member * loopMember, NewColInformation * newColInfo){
    assert(dec);
    assert(reducedMember);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));

    assert(reducedMember->numPathArcs != 0);
    bool createPathSeries = reducedMember->numPathArcs > 1;
    bool convertOriginal = reducedMember->numPathArcs == getNumMemberArcs(dec,member) -1;

    //TODO: for now very elaborate to first get logic correct. Can probably merge some branches later...
    if(!createPathSeries && !convertOriginal){
        spqr_member parallel;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &parallel));
        path_arc_id pathArcId = reducedMember->firstPathArc;
        assert(pathArcIsValid(pathArcId));
        spqr_arc marked = newCol->pathArcs[pathArcId].arc;
        assert(marked != markerToParent(dec,member));//TODO: handle this case later
        moveArcToNewMember(dec,marked,member,parallel);
        bool reversed = arcIsReversedNonRigid(dec,marked) ;
        MATREC_CALL(createMarkerPair(dec,member,parallel,true, reversed,reversed));
        *loopMember = parallel;

        if(reversed == reducedMember->pathBackwards){
            setTerminalReversed(newColInfo,!reversed);
        }else{
            setTerminalReversed(newColInfo,reversed);
        }

        setTerminalMember(newColInfo,*loopMember);
        return MATREC_OKAY;
    }
    if(!createPathSeries && convertOriginal){
        //only one path arc; we are in a loop; no need to change anything
        assert(getNumMemberArcs(dec,member) == 2);
        assert(reducedMember->numPathArcs == 1);
        *loopMember = member;
        changeLoopToParallel(dec,member);

        path_arc_id pathArcId = reducedMember->firstPathArc;
        spqr_arc marked = newCol->pathArcs[pathArcId].arc;
        //The 'reversed' field has different meaning for parallels, so we need to change orientation when converting to a parallel
        arcFlipReversed(dec,marked);

        setTerminalReversed(newColInfo,reducedMember->pathBackwards);
        setTerminalMember(newColInfo,*loopMember);
        return MATREC_OKAY;
        //TODO: this seems too simple, are we relying on assumption that all column edges are forward?
    }
    if(createPathSeries && !convertOriginal){
        spqr_member pathMember;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

        path_arc_id pathArcId = reducedMember->firstPathArc;
        bool parentMoved = false;
        while(pathArcIsValid(pathArcId)){
            spqr_arc pathArc = newCol->pathArcs[pathArcId].arc;
            pathArcId = newCol->pathArcs[pathArcId].nextMember;
            if(pathArc == markerToParent(dec,member)){
                parentMoved = true;
            }
            moveArcToNewMember(dec,pathArc,member,pathMember);
        }

        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, loopMember));

        if(!parentMoved){
            MATREC_CALL(createMarkerPair(dec,member,*loopMember,true,false,false));
            MATREC_CALL(createMarkerPair(dec,*loopMember,pathMember,true,false,true));
        }else{
            MATREC_CALL(createMarkerPair(dec,pathMember,*loopMember,false,false,true));
            MATREC_CALL(createMarkerPair(dec,*loopMember,member,false,false,false));
        }

        setTerminalReversed(newColInfo,!reducedMember->pathBackwards);
        setTerminalMember(newColInfo,*loopMember);
        return MATREC_OKAY;
    }
    assert(createPathSeries && convertOriginal);
    //There's one exception in this case
    //if the single unmarked (column) marker is a parent or child marker to a parallel member, we add the edge there
    {
        spqr_member adjacentMember = SPQR_INVALID_MEMBER;
        spqr_arc adjacentMarker = SPQR_INVALID_ARC;
        spqr_arc memberMarker = SPQR_INVALID_ARC;
        spqr_arc firstArc = getFirstMemberArc(dec, reducedMember->member);
        spqr_arc arc = firstArc;
        do {
            if (!newCol->arcInPath[arc]) {
                if (arc == markerToParent(dec, reducedMember->member)) {
                    adjacentMember = findMemberParent(dec, reducedMember->member);
                    adjacentMarker = markerOfParent(dec, reducedMember->member);
                    memberMarker = arc;
                } else if (arcIsMarker(dec, arc)) {
                    adjacentMember = findArcChildMember(dec, arc);
                    adjacentMarker = markerToParent(dec, adjacentMember);
                    memberMarker = arc;
                }

                break; //There is only a singular such arc
            }
            arc = getNextMemberArc(dec, arc);
        } while (arc != firstArc);

        if (SPQRmemberIsValid(adjacentMember)) {
            SPQRMemberType adjacentType = getMemberType(dec, adjacentMember);
            if (adjacentType == SPQR_MEMBERTYPE_PARALLEL) {
                //Figure out if the markers are the same or opposite orientations
                //If they are the same, we can proceed as normal, otherwise, we need to flip the placed edge
                bool markersHaveSameOrientation =
                        arcIsReversedNonRigid(dec, adjacentMarker) == arcIsReversedNonRigid(dec, memberMarker);
                setTerminalReversed(newColInfo, reducedMember->pathBackwards == markersHaveSameOrientation);
                setTerminalMember(newColInfo, adjacentMember);
                return MATREC_OKAY;
            }
        }
    }

    spqr_member pathMember;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

    path_arc_id pathArcId = reducedMember->firstPathArc;
    bool parentMoved = false;
    while(pathArcIsValid(pathArcId)){
        spqr_arc pathArc = newCol->pathArcs[pathArcId].arc;
        pathArcId = newCol->pathArcs[pathArcId].nextMember;
        if(pathArc == markerToParent(dec,member)){
            parentMoved = true;
        }
        moveArcToNewMember(dec,pathArc,member,pathMember);
    }
    if(parentMoved){
        MATREC_CALL(createMarkerPair(dec,pathMember,member,false,false,false));
    }else{
        MATREC_CALL(createMarkerPair(dec,member,pathMember,true,false,false));
    }

    changeLoopToParallel(dec,member);

    *loopMember = member;
    setTerminalReversed(newColInfo,reducedMember->pathBackwards);
    setTerminalMember(newColInfo,*loopMember);
    return MATREC_OKAY;
}


static MATREC_ERROR splitSeriesMerging(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                                     MATRECColReducedMember * reducedMember,
                                     spqr_member member,
                                     spqr_arc * pathRepresentative,
                                     spqr_arc * nonPathRepresentative,
                                     spqr_arc exceptionArc1,
                                     spqr_arc exceptionArc2){
    assert(dec);
    assert(reducedMember);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec,member));

    int numExceptionArcs = (exceptionArc1 == SPQR_INVALID_ARC ? 0 : 1) + (exceptionArc2 == SPQR_INVALID_ARC ? 0 : 1);
    int numNonPathArcs = getNumMemberArcs(dec,member) - reducedMember->numPathArcs - numExceptionArcs;
    bool createPathSeries = reducedMember->numPathArcs > 1;
    //If this holds, there are 2 or more non-parent marker non-path arcs
    bool createNonPathSeries = numNonPathArcs > 1;
    assert(exceptionArc1 == SPQR_INVALID_ARC || !newCol->arcInPath[exceptionArc1]);
    assert(exceptionArc2 == SPQR_INVALID_ARC || !newCol->arcInPath[exceptionArc2]);

    if(createPathSeries){
        spqr_member pathMember;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

        path_arc_id pathArcId = reducedMember->firstPathArc;
        bool parentMoved = false;
        while(pathArcIsValid(pathArcId)){
            spqr_arc pathArc = newCol->pathArcs[pathArcId].arc;
            pathArcId = newCol->pathArcs[pathArcId].nextMember;
            assert(pathArc != exceptionArc1 && pathArc != exceptionArc2);
            parentMoved = parentMoved || markerToParent(dec,member) == pathArc;
            moveArcToNewMember(dec,pathArc,member,pathMember);
        }
        assert(getNumMemberArcs(dec,pathMember) >= 2);

        spqr_arc ignored;
        bool inOldReversed = true;
        bool inNewReversed = false;
        if(parentMoved){
            MATREC_CALL(createMarkerPairWithReferences(dec,pathMember,member,false,inNewReversed,inOldReversed,
                                                     &ignored,pathRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,pathMember,true,inOldReversed,inNewReversed,
                                                     pathRepresentative,&ignored));
        }
    }else{
        if(pathArcIsValid(reducedMember->firstPathArc)){
            *pathRepresentative = newCol->pathArcs[reducedMember->firstPathArc].arc;
        }else{
            *pathRepresentative = SPQR_INVALID_ARC;
        }
    }

    if(createNonPathSeries){
        spqr_member nonPathMember;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &nonPathMember));

        spqr_arc arc = getFirstMemberArc(dec, member);
        bool parentMoved = false;
        bool canStop = false; //hack when the first arc is moved in the below loop to prevent that we immediately terminate
        do{
            spqr_arc nextArc = getNextMemberArc(dec, arc);
            if(arc != *pathRepresentative && arc != exceptionArc1 && arc != exceptionArc2){
                parentMoved = parentMoved || markerToParent(dec,member) == arc;
                moveArcToNewMember(dec,arc,member,nonPathMember);
            }else{
                canStop = true;
            }
            arc = nextArc;
            if(canStop && arc == getFirstMemberArc(dec,member)){
                break;
            }
        }while(true);
        assert(getNumMemberArcs(dec,nonPathMember) >= 2);
        bool representativeIsTree = !arcIsTree(dec,exceptionArc1);
        if(SPQRarcIsValid(exceptionArc2)){
            representativeIsTree = representativeIsTree || !arcIsTree(dec,exceptionArc2);
        }
        spqr_arc ignored;
        bool inOldReversed = true;
        bool inNewReversed = false;
        if(parentMoved){
            MATREC_CALL(createMarkerPairWithReferences(dec,nonPathMember,member,!representativeIsTree,
                                                     inNewReversed,inOldReversed,&ignored,nonPathRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,nonPathMember,representativeIsTree,
                                                     inOldReversed,inNewReversed,nonPathRepresentative,&ignored));
        }
    }else{
        *nonPathRepresentative = SPQR_INVALID_ARC;
        if(numNonPathArcs != 0) {
            spqr_arc firstArc = getFirstMemberArc(dec, member);
            spqr_arc arc = firstArc;
            do {
                if (arc != *pathRepresentative && arc != exceptionArc1 && arc != exceptionArc2) {
                    *nonPathRepresentative = arc;
                    break;
                }
                arc = getNextMemberArc(dec, arc);
            } while (arc != firstArc);
            assert(*nonPathRepresentative != SPQR_INVALID_ARC);
        }
    }

    return MATREC_OKAY;
}

static MATREC_ERROR transformFirstPathMember(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                                           reduced_member_id reducedMember, NewColInformation * newColInfo,
                                           spqr_arc * representativeArc,
                                           spqr_member * mergedMember ){
    spqr_member member = newCol->reducedMembers[reducedMember].member;
    SPQRMemberType type = getMemberType(dec,member);
    if(type == SPQR_MEMBERTYPE_RIGID){
        //The nodes are already created, we only need to assign the correct start/end node
        switch(newCol->reducedMembers[reducedMember].pathType){
            case INTO_HEAD:
            case INTO_TAIL:
                setTerminalTail(newColInfo,newCol->reducedMembers[reducedMember].rigidPathStart);
                break;
            case OUT_HEAD:
            case OUT_TAIL:
                setTerminalHead(newColInfo,newCol->reducedMembers[reducedMember].rigidPathEnd);
                break;
        }
        *representativeArc = findArcSign(dec,newCol->reducedMembers[reducedMember].pathTargetArc).representative;
        *mergedMember = member;

        return MATREC_OKAY;
    }
    assert(type == SPQR_MEMBERTYPE_SERIES);
    //Split off sets of multiple path non-path edges so that the series has exactly 3 edges

    spqr_arc target = newCol->reducedMembers[reducedMember].pathTargetArc;
    MATRECColReducedMember * redMem = &newCol->reducedMembers[reducedMember];
    spqr_arc pathRepresentative = SPQR_INVALID_ARC;
    spqr_arc nonPathRepresentative = SPQR_INVALID_ARC;
    MATREC_CALL(splitSeriesMerging(dec,newCol,redMem,member,&pathRepresentative,&nonPathRepresentative,target,SPQR_INVALID_ARC));

    assert(target != pathRepresentative
    && target != nonPathRepresentative && pathRepresentative != nonPathRepresentative);
    assert(SPQRarcIsValid(pathRepresentative) && SPQRarcIsValid(nonPathRepresentative) && SPQRarcIsValid(target));
    assert(getNumMemberArcs(dec,member) == 3);

    //Create nodes
    spqr_node a = SPQR_INVALID_NODE;
    spqr_node b = SPQR_INVALID_NODE;
    spqr_node c = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec, &a));
    MATREC_CALL(createNode(dec, &b));
    MATREC_CALL(createNode(dec, &c));

    // a -- b
    //  \  /
    //   c

    //Set arc nodes
    //Set target from b to c,
    bool targetReversed = arcIsReversedNonRigid(dec,target);
    setArcHeadAndTail(dec, target, c, b);

    MemberPathType pathType = newCol->reducedMembers[reducedMember].pathType;
    assert(pathType == INTO_HEAD || pathType == OUT_HEAD);
    if(arcIsReversedNonRigid(dec,pathRepresentative) == targetReversed){
        setArcHeadAndTail(dec,pathRepresentative,a,c);
    }else{
        setArcHeadAndTail(dec,pathRepresentative,c,a);
    }
    if(arcIsReversedNonRigid(dec,nonPathRepresentative) == targetReversed){
        setArcHeadAndTail(dec,nonPathRepresentative,b,a);
    }else{
        setArcHeadAndTail(dec,nonPathRepresentative,a,b);
    }
    //setup signed union find; all arcs are placed are not reversed. We pick an arbitrary arc as 'root' arc for this skeleton
    arcSetReversed(dec,target,false);
    arcSetReversed(dec,pathRepresentative,false);
    arcSetReversed(dec,nonPathRepresentative,false);
    arcSetRepresentative(dec,target,SPQR_INVALID_ARC);
    arcSetRepresentative(dec,pathRepresentative,target);
    arcSetRepresentative(dec,nonPathRepresentative,target);
    *representativeArc = target;

    if(pathType == INTO_HEAD){
        setTerminalTail(newColInfo,a);
    }else{
        setTerminalHead(newColInfo,a);
    }

    *mergedMember = member;

    return MATREC_OKAY;
}
static MATREC_ERROR transformAndMergeParallel(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                                            reduced_member_id current, reduced_member_id next,
                                            spqr_member nextMember, bool nextIsParent,
                                            spqr_arc * representativeArc,
                                            spqr_member * mergedMember){
    //Split off edges not in current subtree to form one parallel (or one edge)
    spqr_member childParallel = INVALID_REDUCED_MEMBER;
    spqr_arc source = newCol->reducedMembers[next].pathSourceArc;
    spqr_arc target = newCol->reducedMembers[next].pathTargetArc;
    if(getNumMemberArcs(dec,nextMember) > 3){
        MATREC_CALL(splitParallel(dec, nextMember, source, target, &childParallel));
        nextMember = childParallel;
        newCol->reducedMembers[next].member = childParallel;
    }
    assert(getNumMemberArcs(dec,nextMember) == 3);

    //TODO: we do some extra work by creating nodes

    spqr_node sourceHead = SPQR_INVALID_NODE;
    spqr_node sourceTail = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec,&sourceHead));
    MATREC_CALL(createNode(dec,&sourceTail));

    //set edge nodes and arc union-find data
    {
        spqr_arc firstArc = getFirstMemberArc(dec,nextMember);
        spqr_arc arc = firstArc;

        bool sourceReversed = arcIsReversedNonRigid(dec,source);
        do{
            if(arcIsReversedNonRigid(dec,arc) == sourceReversed){
                setArcHeadAndTail(dec,arc,sourceHead,sourceTail);
            }else{
                setArcHeadAndTail(dec,arc,sourceTail,sourceHead);
            }
            arcSetRepresentative(dec,arc,source);
            arcSetReversed(dec,arc,false);

            arc = getNextMemberArc(dec,arc);
        }while(arc != firstArc);

        arcSetRepresentative(dec,source,SPQR_INVALID_ARC);
    }

    //fix arc orientations of members; we cannot reflect for parallels.
    *representativeArc = mergeArcSigns(dec,*representativeArc,source,false);

    spqr_member newMergedMember = SPQR_INVALID_MEMBER;
    if(nextIsParent){
        mergeGivenMemberIntoParent(dec,*mergedMember,nextMember,
                                   source,newCol->reducedMembers[current].pathTargetArc
                                   ,true,&newMergedMember);
    }else{
        mergeGivenMemberIntoParent(dec,nextMember,*mergedMember,
                                   newCol->reducedMembers[current].pathTargetArc,source ,true,&newMergedMember);
    }
    * mergedMember = newMergedMember;

    return MATREC_OKAY;
}
static MATREC_ERROR transformAndMergeSeries(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                                          reduced_member_id current, reduced_member_id next,
                                          spqr_member nextMember, bool nextIsParent,
                                          spqr_arc * representativeArc,
                                          spqr_member * mergedMember,
                                          NewColInformation * info){

    MATRECColReducedMember * redMem = &newCol->reducedMembers[next];
    spqr_arc source = redMem->pathSourceArc;
    spqr_arc target = redMem->pathTargetArc;

    spqr_arc pathRepresentative = SPQR_INVALID_ARC;
    spqr_arc nonPathRepresentative = SPQR_INVALID_ARC;
    MATREC_CALL(splitSeriesMerging(dec,newCol,redMem,nextMember,&pathRepresentative,&nonPathRepresentative,source,target));
    //After splitting there is the following possibilities for nodes a-d:
    //(a)-source-(b)-path-(c)-target-(d)-nonpath-(a)
    //(a)-source-(b)-path-(c)-target-(d==a)
    //(a)-source-(b)=(c)-target-(d)-nonpath-(a)
    //(a)-source-(b)-path-(c)=(d) -nonpath-(a)
    //Note that the given arc is always between the same nodes
    assert(getNumMemberArcs(dec,nextMember) == 3 || getNumMemberArcs(dec,nextMember) == 4);
    assert(pathRepresentative != source && nonPathRepresentative != source && (SPQRarcIsInvalid(target) || (
            target != pathRepresentative && target != nonPathRepresentative
            )));
    spqr_node a = SPQR_INVALID_NODE;
    spqr_node b = SPQR_INVALID_NODE;
    spqr_node c = SPQR_INVALID_NODE;
    spqr_node d = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec, &a));
    MATREC_CALL(createNode(dec, &b));
    if (SPQRarcIsValid(pathRepresentative)) {
        MATREC_CALL(createNode(dec, &c));
    } else {
        c = b;
    }
    bool hasNonPath = SPQRarcIsValid(nonPathRepresentative);
    bool hasTarget = SPQRarcIsValid(target);
    if (hasNonPath && hasTarget) {
        MATREC_CALL(createNode(dec, &d));
    } else {
        if(hasNonPath){
            d = c;
        }else{
            d = a;
        }
    }

    bool sourceReversed = arcIsReversedNonRigid(dec,source);
    bool pathStartInHead = isHead(newCol->reducedMembers[current].pathType);
    if(pathStartInHead){
        setArcHeadAndTail(dec,source,b,a);
    }else{
        setArcHeadAndTail(dec,source,a,b);
    }
    if(SPQRarcIsValid(pathRepresentative)){
        if((arcIsReversedNonRigid(dec,pathRepresentative) == sourceReversed) == pathStartInHead){
            setArcHeadAndTail(dec,pathRepresentative,c,b);
        }else{
            setArcHeadAndTail(dec,pathRepresentative,b,c);
        }
        arcSetReversed(dec,pathRepresentative,false);
        arcSetRepresentative(dec,pathRepresentative,source);
    }
    if(hasTarget){
        if((arcIsReversedNonRigid(dec,target) == sourceReversed) == pathStartInHead){
            setArcHeadAndTail(dec,target,d,c);
        }else{
            setArcHeadAndTail(dec,target,c,d);
        }
        arcSetReversed(dec,target,false);
        arcSetRepresentative(dec,target,source);
    }
    if(hasNonPath){
        if((arcIsReversedNonRigid(dec,nonPathRepresentative) == sourceReversed) == pathStartInHead){
            setArcHeadAndTail(dec,nonPathRepresentative,a,d);
        }else{
            setArcHeadAndTail(dec,nonPathRepresentative,d,a);
        }
        arcSetReversed(dec,nonPathRepresentative,false);
        arcSetRepresentative(dec,nonPathRepresentative,source);
    }

    arcSetReversed(dec,source,false);
    arcSetRepresentative(dec,source,SPQR_INVALID_ARC);

    //fix arc orientations of members; we cannot reflect for series

    spqr_member newMergedMember = SPQR_INVALID_MEMBER;
    if(nextIsParent){
        mergeGivenMemberIntoParent(dec,*mergedMember,nextMember,
                                   source,newCol->reducedMembers[current].pathTargetArc
                ,true,&newMergedMember);
    }else{
        mergeGivenMemberIntoParent(dec,nextMember,*mergedMember,
                                   newCol->reducedMembers[current].pathTargetArc,source ,true,&newMergedMember);
    }
    * mergedMember = newMergedMember;

    *representativeArc = mergeArcSigns(dec,*representativeArc,source,false);
    if(!hasTarget){
        //We are in the last node; finish the path
        setTerminalReversed(info,false);
        if(isInto(newCol->reducedMembers[current].pathType)) {
            setTerminalHead(info, c);
        }else{
            setTerminalTail(info, c);
        }
        setTerminalMember(info,*mergedMember);
        setTerminalRepresentative(info,*representativeArc);
    }
    return MATREC_OKAY;
}
static MATREC_ERROR transformAndMergeRigid(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                                         reduced_member_id current, reduced_member_id next,
                                         spqr_member nextMember, bool nextIsParent,
                                         spqr_arc * representativeArc,
                                         spqr_member * mergedMember,
                                         NewColInformation * info){
    MATRECColReducedMember * redMem = &newCol->reducedMembers[next];
    spqr_arc source = redMem->pathSourceArc;
    spqr_arc sourceRepresentative = findArcSign(dec,source).representative;


    spqr_member newMergedMember = SPQR_INVALID_MEMBER;

    if(nextIsParent){
        //TODO: head-to-head here, or is reversal above sufficient?
        mergeGivenMemberIntoParent(dec,*mergedMember,nextMember,
                                   source,newCol->reducedMembers[current].pathTargetArc
                ,!redMem->reverseArcs,&newMergedMember);
    }else{
        mergeGivenMemberIntoParent(dec,nextMember,*mergedMember,
                                   newCol->reducedMembers[current].pathTargetArc,source ,!redMem->reverseArcs,&newMergedMember);
    }

    * mergedMember = newMergedMember;

    *representativeArc = mergeArcSigns(dec,*representativeArc,sourceRepresentative,redMem->reverseArcs);

    if(SPQRarcIsInvalid(redMem->pathTargetArc)){
        //We are in the last node; finish the path
        setTerminalReversed(info,false);
        if(isInto(newCol->reducedMembers[current].pathType)) {
            if(redMem->reverseArcs){
                setTerminalHead(info,redMem->rigidPathStart);
            }else{
                setTerminalHead(info, redMem->rigidPathEnd);
            }
        }else{
            if(redMem->reverseArcs){
                setTerminalTail(info,redMem->rigidPathEnd);
            }else{
                setTerminalTail(info, redMem->rigidPathStart);
            }
        }
        setTerminalMember(info,*mergedMember);
        setTerminalRepresentative(info,*representativeArc);
    }

    return MATREC_OKAY;
}
static MATREC_ERROR transformAndMerge(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                                           reduced_member_id current, reduced_member_id next,
                                           spqr_arc * representativeArc,
                                           spqr_member * mergedMember,
                                           bool nextIsParent, NewColInformation * info){
    spqr_member nextMember = newCol->reducedMembers[next].member;
    switch(getMemberType(dec,nextMember)){
        case SPQR_MEMBERTYPE_RIGID:{
            MATREC_CALL(transformAndMergeRigid(dec,newCol,current,next,nextMember,nextIsParent,
                                             representativeArc,mergedMember,info));
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:{
            MATREC_CALL(transformAndMergeParallel(dec,newCol,current,next,nextMember,nextIsParent,
                                                representativeArc, mergedMember));
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:{
            MATREC_CALL(transformAndMergeSeries(dec,newCol,current,next,nextMember,nextIsParent,
                                                representativeArc, mergedMember,info));
            break;
        }
        case SPQR_MEMBERTYPE_LOOP:
        case SPQR_MEMBERTYPE_UNASSIGNED:
        {
            assert(false);
            break;
        }
    }
    return MATREC_OKAY;
}

static MATREC_ERROR transformPath(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition * newCol,
                                MATRECColReducedComponent * component,
                                NewColInformation * newColInfo){
    //Realize first member
    reduced_member_id firstPathMember = component->pathEndMembers[0];

    spqr_arc representativeArc = SPQR_INVALID_ARC;
    spqr_member mergedMember = SPQR_INVALID_MEMBER;
    MATREC_CALL(transformFirstPathMember(dec,newCol,firstPathMember,newColInfo,&representativeArc,&mergedMember));
    //Iteratively call function which realizes next member and merges them together.
    reduced_member_id current = firstPathMember;
    reduced_member_id next = newCol->reducedMembers[current].nextPathMember;
    bool nextIsParent = newCol->reducedMembers[current].nextPathMemberIsParent;
    while(reducedMemberIsValid(next)){
        MATREC_CALL(transformAndMerge(dec,newCol,current,next,&representativeArc,&mergedMember,nextIsParent,newColInfo));
        current = next;
        next = newCol->reducedMembers[current].nextPathMember;
        nextIsParent = newCol->reducedMembers[current].nextPathMemberIsParent;
    }
    return MATREC_OKAY;
}

static MATREC_ERROR columnTransformSingleParallel(MATRECNetworkColumnAddition *newCol,
                                                reduced_member_id reducedMemberId, spqr_member member,
                                                NewColInformation *newColInfo){
    MATRECColReducedMember * reducedMember = &newCol->reducedMembers[reducedMemberId];
    assert(pathArcIsValid(reducedMember->firstPathArc) && reducedMember->numPathArcs == 1);
    //The new arc can be placed in parallel; just add it to this member
    setTerminalReversed(newColInfo,reducedMember->pathBackwards);
    setTerminalMember(newColInfo,member);
    return MATREC_OKAY;
}
static MATREC_ERROR columnTransformSingleSeries(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol,
                                                reduced_member_id reducedMemberId, spqr_member member,
                                                NewColInformation *newColInfo){

    if(getNumMemberArcs(dec,member) == 1){
        newColInfo->member = member;
        newColInfo->reversed = newCol->arcInPathReversed[getFirstMemberArc(dec,member)];
        return MATREC_OKAY;
    }
    //Isolated single cycle
    spqr_member loopMember;
    MATRECColReducedMember *reducedMember = &newCol->reducedMembers[reducedMemberId];
    MATREC_CALL(splitSeries(dec,newCol,reducedMember,member,&loopMember,newColInfo));

    return MATREC_OKAY;
}
static MATREC_ERROR columnTransformSingleRigid(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol,
                                             reduced_member_id reducedMemberId, spqr_member member,
                                             NewColInformation *newColInfo){
    assert(dec);
    assert(newCol);
    assert(newColInfo);
    assert(reducedMemberIsValid(reducedMemberId));
    //The path is already computed, so we can simply take the start and end nodes.
    //However, there is one exception, which is that an arc connecting these two nodes already exists in the member
    //If so, we create a new parallel member with the new arc and this member, unless the existing arc already points
    //to a parallel member
    MATRECColReducedMember * reducedMember = &newCol->reducedMembers[reducedMemberId];
    assert(SPQRnodeIsValid(reducedMember->rigidPathStart) && SPQRnodeIsValid(reducedMember->rigidPathEnd));
    {
        //TODO: this search is inefficient, technically
        spqr_arc existingArcWithPath = SPQR_INVALID_ARC;
        spqr_arc firstArc = getFirstNodeArc(dec,reducedMember->rigidPathStart);
        spqr_arc arc = firstArc;
        bool pathInSameDirection = false;
        do{
            spqr_node head = findArcHead(dec,arc);
            spqr_node tail = findArcTail(dec,arc);
            spqr_node other = head == reducedMember->rigidPathStart ? tail : head;
            if(other == reducedMember->rigidPathEnd){
                existingArcWithPath = arc;
                pathInSameDirection = (head == other) != findArcSign(dec,existingArcWithPath).reversed;
                break;
            }
            arc = getNextNodeArc(dec,arc,reducedMember->rigidPathStart);
        }while(arc != firstArc);
        if(SPQRarcIsValid(existingArcWithPath)){
            bool isParent = false;
            spqr_member adjacentMember = arcIsMarker(dec,existingArcWithPath) ? findArcChildMember(dec,existingArcWithPath) : SPQR_INVALID_MEMBER;
            if(existingArcWithPath == markerToParent(dec,member)){
                adjacentMember = findMemberParent(dec,member);
                isParent = true;
            }
            if(SPQRmemberIsValid(adjacentMember) && getMemberType(dec,adjacentMember) == SPQR_MEMBERTYPE_PARALLEL){
                spqr_arc parallelMarker = isParent ? markerOfParent(dec,member) : markerToParent(dec,adjacentMember);
                bool markerReversed = arcIsReversedNonRigid(dec,parallelMarker);
                setTerminalMember(newColInfo,adjacentMember);
                setTerminalReversed(newColInfo,markerReversed == pathInSameDirection);
            }else{
                //create a new parallel and move the edge there
                //This is a bit painful, because we cannot actually remove edges because of the union-find data structure
                //So what we do instead, is convert the current edge to a marker edge, and 'duplicate'
                //it in the new parallel member, and add the new marker there too, manually
                spqr_member adjacentParallel = SPQR_INVALID_MEMBER;
                MATREC_CALL(createMember(dec,SPQR_MEMBERTYPE_PARALLEL,&adjacentParallel));
                //'duplicate' a new arc in the parallel to be the current arc
                spqr_arc duplicate = SPQR_INVALID_ARC;
                spqr_element element = arcGetElement(dec,existingArcWithPath);
                if( element != MARKER_COLUMN_ELEMENT && element != MARKER_ROW_ELEMENT){
                    if(SPQRelementIsColumn(element)){
                        MATREC_CALL(createColumnArc(dec,adjacentParallel,&duplicate, SPQRelementToColumn(element),false));
                    }else{
                        MATREC_CALL(createRowArc(dec,adjacentParallel,&duplicate, SPQRelementToRow(element),false));
                    }
                }else if(isParent){
                    MATREC_CALL(createParentMarker(dec,adjacentParallel, arcIsTree(dec,existingArcWithPath),adjacentMember,
                                                 markerOfParent(dec,member)
                                                 ,&duplicate,false));
                }else{
                    MATREC_CALL(createChildMarker(dec,adjacentParallel,adjacentMember,arcIsTree(dec,existingArcWithPath),&duplicate,false));
                    dec->members[adjacentMember].parentMember = adjacentParallel;
                    dec->members[adjacentMember].markerOfParent = duplicate;
                }
                //Create the other marker edge
                spqr_arc parallelMarker = SPQR_INVALID_ARC;
                if(isParent){
                    MATREC_CALL(createChildMarker(dec,adjacentParallel,member,!arcIsTree(dec,existingArcWithPath),
                                                &parallelMarker,false));
                }else{
                    MATREC_CALL(createParentMarker(dec,adjacentParallel,!arcIsTree(dec,existingArcWithPath),
                                                 member,existingArcWithPath,&parallelMarker,false));
                }

                //Change the existing edge to a marker
                if(isParent){
                    assert(markerToParent(dec,member) == existingArcWithPath);
                    dec->arcs[markerOfParent(dec,member)].childMember = adjacentParallel;
                    dec->members[member].parentMember = adjacentParallel;
                    dec->members[member].markerToParent = existingArcWithPath;
                    dec->members[member].markerOfParent = parallelMarker;
                    dec->arcs[existingArcWithPath].element =  arcIsTree(dec,existingArcWithPath) ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;;
                    dec->arcs[existingArcWithPath].childMember = adjacentParallel;

                }else{
                    dec->arcs[existingArcWithPath].element = arcIsTree(dec,existingArcWithPath) ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;
                    dec->arcs[existingArcWithPath].childMember = adjacentParallel;
                }

                setTerminalMember(newColInfo,adjacentParallel);
                setTerminalReversed(newColInfo,!pathInSameDirection);
            }
            return MATREC_OKAY;
        }
    }

    setTerminalMember(newColInfo,member);
    setTerminalReversed(newColInfo,false);
    setTerminalTail(newColInfo,reducedMember->rigidPathStart);
    setTerminalHead(newColInfo,reducedMember->rigidPathEnd);
    setTerminalRepresentative(newColInfo,findArcSign(dec,newCol->pathArcs[reducedMember->firstPathArc].arc).representative);
    return MATREC_OKAY;
}
static MATREC_ERROR transformComponent(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol, MATRECColReducedComponent * component, NewColInformation * newColInfo){
    assert(dec);
    assert(newCol);
    assert(component);
    assert(newColInfo);

    if(newCol->reducedMembers[component->root].numChildren == newCol->reducedMembers[component->root].numPropagatedChildren){
        //No merging necessary, only a single component
        reduced_member_id reducedMember = component->root;
        assert(reducedMemberIsValid(reducedMember));
        spqr_member member = newCol->reducedMembers[reducedMember].member;
        SPQRMemberType type = getMemberType(dec, member);

        switch(type) {
            case SPQR_MEMBERTYPE_RIGID: {
                MATREC_CALL(columnTransformSingleRigid(dec,newCol,reducedMember,member,newColInfo));
                break;
            }
            case SPQR_MEMBERTYPE_PARALLEL: {
                MATREC_CALL(columnTransformSingleParallel( newCol, reducedMember, member, newColInfo));
                break;
            }
            case SPQR_MEMBERTYPE_LOOP:
            case SPQR_MEMBERTYPE_SERIES: {
                MATREC_CALL(columnTransformSingleSeries(dec, newCol, reducedMember, member, newColInfo));
                break;
            }
            default: {
                assert(false);
                break;
            }
        }
        return MATREC_OKAY;
    }
    // Otherwise, the reduced members form a path which can be merged into a single component of type R
    MATREC_CALL(transformPath(dec,newCol,component,newColInfo));

    return MATREC_OKAY;
}

MATREC_ERROR MATRECNetworkColumnAdditionAdd(MATRECNetworkDecomposition *dec, MATRECNetworkColumnAddition *newCol){
    assert(dec);
    assert(newCol);

    if(newCol->numReducedComponents == 0){
        spqr_member member;
        MATREC_CALL(createStandaloneSeries(dec,newCol->newRowArcs,newCol->newRowArcReversed,
                                         newCol->numNewRowArcs,newCol->newColIndex,&member));
    }else if(newCol->numReducedComponents == 1){
        NewColInformation information = emptyNewColInformation();
        MATREC_CALL(transformComponent(dec,newCol,&newCol->reducedComponents[0],&information));
        assert(memberIsRepresentative(dec,information.member));
        if(newCol->numNewRowArcs == 0){
            spqr_arc colArc = SPQR_INVALID_ARC;
            MATREC_CALL(createColumnArc(dec,information.member,&colArc,newCol->newColIndex,information.reversed));
            if(SPQRnodeIsValid(information.head)){
                assert(SPQRnodeIsValid(information.tail));
                assert(SPQRarcIsValid(information.representative));
                //TODO: check if finds are necessary
                setArcHeadAndTail(dec,colArc,findNode(dec,information.head),findNode(dec,information.tail));
                arcSetRepresentative(dec,colArc,information.representative);
                arcSetReversed(dec,colArc,information.reversed != arcIsReversedNonRigid(dec,information.representative));
            }
        }else{
            spqr_member newSeries = SPQR_INVALID_MEMBER;
            MATREC_CALL(createConnectedSeries(dec,newCol->newRowArcs,newCol->newRowArcReversed,newCol->numNewRowArcs,newCol->newColIndex,&newSeries));
            spqr_arc markerArc = SPQR_INVALID_ARC;
            spqr_arc ignore = SPQR_INVALID_ARC;
            MATREC_CALL(createMarkerPairWithReferences(dec,information.member,newSeries,false,information.reversed,true,&markerArc,&ignore));
            if(SPQRnodeIsValid(information.head)){
                assert(SPQRnodeIsValid(information.tail));
                assert(SPQRarcIsValid(information.representative));
                //TODO: check if finds are necessary
                setArcHeadAndTail(dec,markerArc,findNode(dec,information.head),findNode(dec,information.tail));
                arcSetRepresentative(dec,markerArc,information.representative);
                arcSetReversed(dec,markerArc,information.reversed != arcIsReversedNonRigid(dec,information.representative));
            }
        }
        if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
            assert(getNumMemberArcs(dec,information.member) == 2 || getNumMemberArcs(dec,information.member) == 3);
            if(getNumMemberArcs(dec,information.member) == 3){
                changeLoopToParallel(dec,information.member);
            }
        }
    }else{
#ifndef NDEBUG
        int numDecComponentsBefore = numConnectedComponents(dec);
#endif
        spqr_member newSeries = SPQR_INVALID_MEMBER;
        MATREC_CALL(createConnectedSeries(dec,newCol->newRowArcs,newCol->newRowArcReversed,
                                        newCol->numNewRowArcs,newCol->newColIndex,&newSeries));
        for (int i = 0; i < newCol->numReducedComponents; ++i) {
            NewColInformation information = emptyNewColInformation();
            MATREC_CALL(transformComponent(dec,newCol,&newCol->reducedComponents[i],&information));
            if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
                assert(getNumMemberArcs(dec,information.member) == 1);
                spqr_arc arc = getFirstMemberArc(dec,information.member);
                assert(newCol->arcInPath[arc]);
                moveArcToNewMember(dec,arc,information.member,newSeries);
                arcSetReversed(dec,arc, !newCol->arcInPathReversed[arc]);
                dec->members[information.member].type = SPQR_MEMBERTYPE_UNASSIGNED;
            }else {
                reorderComponent(dec,
                                 information.member); //reorder the subtree so that the newly series member is a parent
                spqr_arc markerArc = SPQR_INVALID_ARC;
                spqr_arc ignore = SPQR_INVALID_ARC;
                MATREC_CALL(createMarkerPairWithReferences(dec, newSeries, information.member, true, information.reversed,
                                                         true, &ignore, &markerArc));
                if (SPQRnodeIsValid(information.head)) {
                    assert(SPQRnodeIsValid(information.tail));
                    assert(SPQRarcIsValid(information.representative));
                    //TODO: check if finds are necessary
                    setArcHeadAndTail(dec, markerArc, findNode(dec, information.head), findNode(dec, information.tail));
                    arcSetRepresentative(dec, markerArc, information.representative);
                    //TODO: is this correct?
                    arcSetReversed(dec, markerArc,
                                   information.reversed == arcIsReversedNonRigid(dec, information.representative));
                }
            }
        }
        decreaseNumConnectedComponents(dec,newCol->numReducedComponents-1);
        assert(numConnectedComponents(dec) == (numDecComponentsBefore - newCol->numReducedComponents + 1));
    }
//    decompositionToDot(stdout,dec,true);
    return MATREC_OKAY;
}

bool MATRECNetworkColumnAdditionRemainsNetwork(MATRECNetworkColumnAddition *newCol){
    return newCol->remainsNetwork;
}


static int min(int a, int b){
    return a < b ? a : b;
}

typedef int cut_arc_id;
#define INVALID_CUT_ARC (-1)

static bool cutArcIsInvalid(const cut_arc_id arc){
    return arc < 0;
}
static bool cutArcIsValid(const cut_arc_id arc){
    return !cutArcIsInvalid(arc);
}

typedef struct { //TODO:test if overhead of pointers is worth it?
    spqr_arc arc;
    spqr_node arcHead;
    spqr_node arcTail;
    cut_arc_id nextMember;
    cut_arc_id nextOverall;
    bool arcReversed;
} CutArcListNode;

typedef enum{
    TYPE_UNDETERMINED = 0,
    TYPE_PROPAGATED = 1,
    TYPE_MERGED = 2,
    TYPE_NOT_NETWORK = 3
} RowReducedMemberType;


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
    spqr_arc nodeArc;
} DFSCallData;

typedef struct{
    children_idx currentChild;
    reduced_member_id id;
} MergeTreeCallData;

typedef struct{
    spqr_node node;
    spqr_arc arc;
} ColorDFSCallData;

typedef struct{
    spqr_arc arc;
    spqr_node node;
    spqr_node parent;
    bool isAP;
} ArticulationPointCallStack;

typedef enum{
    UNCOLORED = 0,
    COLOR_SOURCE = 1,
    COLOR_SINK = 2
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

    cut_arc_id firstCutArc;
    int numCutArcs;

    //For non-rigid members
    spqr_arc splitArc;
    bool splitHead; //Otherwise the tail of this arc is split
    bool otherIsSource; //Otherwise the other end node is part of the sink partition.
    // For non-rigid members this refers to splitArc, for rigid members it refers to articulation arc

    //For rigid members
    spqr_node otherNode;
    spqr_node splitNode;
    bool allHaveCommonNode;
    bool otherNodeSplit;
    bool willBeReversed;
    spqr_arc articulationArc;
    spqr_node coloredNode; //points to a colored node so that we can efficiently zero out the colors again.

} MATRECRowReducedMember;

typedef struct {
    int rootDepth;
    reduced_member_id root;
} MATRECRowReducedComponent;

struct MATRECNetworkRowAdditionImpl {
    bool remainsNetwork;

    MATRECRowReducedMember *reducedMembers;
    int memReducedMembers;
    int numReducedMembers;

    MATRECRowReducedComponent *reducedComponents;
    int memReducedComponents;
    int numReducedComponents;

    MemberInfo *memberInformation;
    int memMemberInformation;
    int numMemberInformation;

    reduced_member_id *childrenStorage;
    int memChildrenStorage;
    int numChildrenStorage;

    CutArcListNode *cutArcs;
    int memCutArcs;
    int numCutArcs;
    cut_arc_id firstOverallCutArc;

    MATREC_row newRowIndex;

    MATREC_col *newColumnArcs;
    bool *newColumnReversed;
    int memColumnArcs;
    int numColumnArcs;

    reduced_member_id *leafMembers;
    int numLeafMembers;
    int memLeafMembers;

    spqr_arc *decompositionColumnArcs;
    bool *decompositionColumnArcReversed;
    int memDecompositionColumnArcs;
    int numDecompositionColumnArcs;

    bool *isArcCut;
    bool *isArcCutReversed;
    int numIsArcCut;
    int memIsArcCut;

    COLOR_STATUS *nodeColors;
    int memNodeColors;

    spqr_node *articulationNodes;
    int numArticulationNodes;
    int memArticulationNodes;

    ArticulationNodeInformation *articulationNodeSearchInfo;
    int memNodeSearchInfo;

    int *crossingPathCount;
    int memCrossingPathCount;

    DFSCallData *intersectionDFSData;
    int memIntersectionDFSData;

    ColorDFSCallData *colorDFSData;
    int memColorDFSData;

    ArticulationPointCallStack *artDFSData;
    int memArtDFSData;

    CreateReducedMembersCallstack *createReducedMembersCallstack;
    int memCreateReducedMembersCallstack;

    int *intersectionPathDepth;
    int memIntersectionPathDepth;

    spqr_node *intersectionPathParent;
    int memIntersectionPathParent;

    MergeTreeCallData *mergeTreeCallData;
    int memMergeTreeCallData;
};

typedef struct {
    spqr_member member;
    spqr_node head;
    spqr_node tail;
    spqr_arc representative;
    bool reversed;
} NewRowInformation;

static NewRowInformation emptyNewRowInformation(void){
    NewRowInformation information;
    information.member = SPQR_INVALID_MEMBER;
    information.head = SPQR_INVALID_NODE;
    information.tail = SPQR_INVALID_NODE;
    information.representative = SPQR_INVALID_ARC;
    information.reversed = false;
    return information;
}

/**
 * Saves the information of the current row and partitions it based on whether or not the given columns are
 * already part of the decomposition.
 */
static MATREC_ERROR newRowUpdateRowInformation(const MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                               const MATREC_row row, const MATREC_col * columns, const double * columnValues,
                                               const size_t numColumns)
{
    newRow->newRowIndex = row;

    newRow->numDecompositionColumnArcs = 0;
    newRow->numColumnArcs = 0;

    for (size_t i = 0; i < numColumns; ++i) {
        spqr_arc columnArc = getDecompositionColumnArc(dec, columns[i]);
        bool reversed = columnValues[i] < 0.0;
        if(SPQRarcIsValid(columnArc)){ //If the arc is the current decomposition: save it in the array
            if(newRow->numDecompositionColumnArcs == newRow->memDecompositionColumnArcs){
                int newNumArcs = newRow->memDecompositionColumnArcs == 0 ? 8 : 2*newRow->memDecompositionColumnArcs; //TODO: make reallocation numbers more consistent with rest?
                newRow->memDecompositionColumnArcs = newNumArcs;
                MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->decompositionColumnArcs,
                                                (size_t) newRow->memDecompositionColumnArcs));
                MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->decompositionColumnArcReversed,
                                                (size_t) newRow->memDecompositionColumnArcs));
            }
            newRow->decompositionColumnArcs[newRow->numDecompositionColumnArcs] = columnArc;
            newRow->decompositionColumnArcReversed[newRow->numDecompositionColumnArcs] = reversed;
            ++newRow->numDecompositionColumnArcs;
        }else{
            //Not in the decomposition: add it to the set of arcs which are newly added with this row.
            if(newRow->numColumnArcs == newRow->memColumnArcs){
                int newNumArcs = newRow->memColumnArcs == 0 ? 8 : 2*newRow->memColumnArcs; //TODO: make reallocation numbers more consistent with rest?
                newRow->memColumnArcs = newNumArcs;
                MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->newColumnArcs,
                                                (size_t)newRow->memColumnArcs));
                MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->newColumnReversed,
                                                (size_t)newRow->memColumnArcs));
            }
            newRow->newColumnArcs[newRow->numColumnArcs] = columns[i];
            newRow->newColumnReversed[newRow->numColumnArcs] = reversed;
            newRow->numColumnArcs++;
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
static reduced_member_id createRowReducedMembersToRoot(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition * newRow, const spqr_member firstMember){
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
            reducedMemberData->numCutArcs = 0;
            reducedMemberData->firstCutArc = INVALID_CUT_ARC;
            reducedMemberData->type = TYPE_UNDETERMINED;
            reducedMemberData->numPropagatedChildren = 0;
            reducedMemberData->articulationArc = SPQR_INVALID_ARC;
            reducedMemberData->splitNode = SPQR_INVALID_NODE;
            reducedMemberData->otherNode = SPQR_INVALID_NODE;
            reducedMemberData->splitArc = SPQR_INVALID_ARC;
            reducedMemberData->splitHead = false;
            reducedMemberData->allHaveCommonNode = false;
            reducedMemberData->otherNodeSplit = false;
            reducedMemberData->willBeReversed = false;
            reducedMemberData->coloredNode = SPQR_INVALID_NODE;

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
 * Construct a smaller sub tree of the decomposition on which the cut arcs lie.
 * @return
 */
static MATREC_ERROR constructRowReducedDecomposition(MATRECNetworkDecomposition* dec, MATRECNetworkRowAddition* newRow){
    //TODO: chop up into more functions
    //TODO: stricter assertions/array bounds checking in this function
#ifndef NDEBUG
    for (int i = 0; i < newRow->memMemberInformation; ++i) {
        assert(reducedMemberIsInvalid(newRow->memberInformation[i].reducedMember));
    }
#endif

    newRow->numReducedComponents = 0;
    newRow->numReducedMembers = 0;
    if(newRow->numDecompositionColumnArcs == 0){ //Early return in case the reduced decomposition will be empty
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
    for (int i = 0; i < newRow->numDecompositionColumnArcs; ++i) {
        assert(i < newRow->memDecompositionColumnArcs);
        spqr_arc arc = newRow->decompositionColumnArcs[i];
        spqr_member arcMember = findArcMember(dec, arc);
        reduced_member_id reducedMember = createRowReducedMembersToRoot(dec,newRow,arcMember);
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
 * Marks an arc as 'cut'. This implies that its cycle in the decomposition must be elongated
 * @param newRow
 * @param arc
 * @param reducedMember
 */
static void createCutArc(MATRECNetworkDecomposition  * dec,
                         MATRECNetworkRowAddition* newRow, const spqr_arc arc, const reduced_member_id reducedMember,
                         bool reversed){
    cut_arc_id cut_arc =  newRow->numCutArcs;

    CutArcListNode * listNode = &newRow->cutArcs[cut_arc];
    listNode->arc = arc;

    listNode->nextMember = newRow->reducedMembers[reducedMember].firstCutArc;
    newRow->reducedMembers[reducedMember].firstCutArc = cut_arc;

    listNode->nextOverall = newRow->firstOverallCutArc;
    newRow->firstOverallCutArc = cut_arc;

    newRow->numCutArcs++;
    newRow->reducedMembers[reducedMember].numCutArcs++;
    assert(newRow->numCutArcs <= newRow->memCutArcs);

    assert(arc < newRow->memIsArcCut);
    newRow->isArcCut[arc] = true;
    newRow->isArcCutReversed[arc] = reversed;

    assert(memberIsRepresentative(dec,newRow->reducedMembers[reducedMember].member));
    if(getMemberType(dec,newRow->reducedMembers[reducedMember].member) == SPQR_MEMBERTYPE_RIGID){

        listNode->arcHead = findEffectiveArcHead(dec,arc);
        listNode->arcTail = findEffectiveArcTail(dec,arc);
        if(reversed){
            swap_ints(&listNode->arcHead,&listNode->arcTail);
        }
        assert(SPQRnodeIsValid(listNode->arcHead) && SPQRnodeIsValid(listNode->arcTail));
    }else{
        listNode->arcHead = SPQR_INVALID_NODE;
        listNode->arcTail = SPQR_INVALID_NODE;
    }

    listNode->arcReversed = reversed;
}

/**
 * Creates all cut arcs within the decomposition for the new row.
 * Note this preallocates memory for cut arcs which may be created by propagation.
 */
static MATREC_ERROR createReducedDecompositionCutArcs(MATRECNetworkDecomposition* dec, MATRECNetworkRowAddition* newRow){
    //Allocate memory for cut arcs
    spqr_arc maxArcID = largestArcID(dec);
    if(maxArcID > newRow->memIsArcCut){
        int newSize = max(maxArcID,2*newRow->memIsArcCut);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->isArcCut,(size_t) newSize));
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->isArcCutReversed,(size_t) newSize));
        for (int i = newRow->memIsArcCut; i < newSize ; ++i) {
            newRow->isArcCut[i] = false;
            newRow->isArcCutReversed[i] = false;
        }
        newRow->memIsArcCut = newSize;
    }
#ifndef NDEBUG
    for (int i = 0; i < newRow->memIsArcCut; ++i) {
        assert(!newRow->isArcCut[i]);
        assert(!newRow->isArcCutReversed[i]);
    }
#endif

    int numNeededArcs = newRow->numDecompositionColumnArcs*4;
    if(numNeededArcs > newRow->memCutArcs){
        int newSize = max(newRow->memCutArcs*2, numNeededArcs);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->cutArcs,(size_t) newSize));
        newRow->memCutArcs = newSize;
    }
    newRow->numCutArcs = 0;
    newRow->firstOverallCutArc = INVALID_CUT_ARC;
    for (int i = 0; i < newRow->numDecompositionColumnArcs; ++i) {
        spqr_arc arc = newRow->decompositionColumnArcs[i];
        spqr_member member = findArcMember(dec, arc);
        reduced_member_id reduced_member = newRow->memberInformation[member].reducedMember;
        assert(reducedMemberIsValid(reduced_member));
        createCutArc(dec,newRow,arc,reduced_member,newRow->decompositionColumnArcReversed[i]);
    }

    return MATREC_OKAY;
}

/**
 * Determines the members of the reduced decomposition which are leafs.
 * This is used in propagation to ensure propagation is only checked for components which have at most one neighbour
 * which is not propagated.
 */
static MATREC_ERROR determineLeafReducedMembers(const MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow) {
    if (newRow->numDecompositionColumnArcs > newRow->memLeafMembers) {
        newRow->memLeafMembers = max(newRow->numDecompositionColumnArcs, 2 * newRow->memLeafMembers);
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
static MATREC_ERROR allocateRigidSearchMemory(const MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow){
    int totalNumNodes = getNumNodes(dec);
    int maxNumNodes = 2*dec->numArcs;
    if(maxNumNodes > newRow->memNodeColors){
        int newSize = max(2*newRow->memNodeColors,maxNumNodes);
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
static void zeroOutColors(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow, const spqr_node firstRemoveNode){
    assert(firstRemoveNode < newRow->memNodeColors);

    newRow->nodeColors[firstRemoveNode] = UNCOLORED;
    ColorDFSCallData * data = newRow->colorDFSData;
    if(newRow->colorDFSData == NULL){
        return;
    }

    data[0].node = firstRemoveNode;
    data[0].arc = getFirstNodeArc(dec,firstRemoveNode);
    if(SPQRarcIsInvalid(data[0].arc)){
        return;
    }

    int depth = 0;

    while(depth >= 0){
        assert(depth < newRow->memColorDFSData);
        ColorDFSCallData * callData = &data[depth];
        spqr_node head = findArcHead(dec, callData->arc);
        spqr_node tail = findArcTail(dec, callData->arc);
        spqr_node otherNode = callData->node == head ? tail : head;
        assert(otherNode < newRow->memNodeColors);
        if(newRow->nodeColors[otherNode] != UNCOLORED){
            callData->arc = getNextNodeArc(dec,callData->arc,callData->node);

            newRow->nodeColors[otherNode] = UNCOLORED;
            ++depth;
            data[depth].node = otherNode;
            data[depth].arc = getFirstNodeArc(dec,otherNode);
            continue;
        }

        callData->arc = getNextNodeArc(dec,callData->arc,callData->node);
        while(depth >= 0 && data[depth].arc == getFirstNodeArc(dec,data[depth].node)){
            --depth;
        }
    }


}
static void cleanUpPreviousIteration(MATRECNetworkDecomposition * dec, MATRECNetworkRowAddition * newRow){
    //zero out coloring information from previous check
    for (int i = 0; i < newRow->numReducedMembers; ++i) {
        if (SPQRnodeIsValid(newRow->reducedMembers[i].coloredNode)) {
            zeroOutColors(dec, newRow, newRow->reducedMembers[i].coloredNode);
        }
    }

#ifndef NDEBUG
    for (int i = 0; i < newRow->memNodeColors; ++i) {
        assert(newRow->nodeColors[i] == UNCOLORED);
    }
#endif

    //For cut arcs: clear them from the array from previous iteration
    cut_arc_id cutArcIdx = newRow->firstOverallCutArc;
    while(cutArcIsValid(cutArcIdx)){
        spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
        cutArcIdx = newRow->cutArcs[cutArcIdx].nextOverall;
        newRow->isArcCut[cutArc] = false;
        newRow->isArcCutReversed[cutArc] = false;
    }
#ifndef NDEBUG
    for (int i = 0; i < newRow->memIsArcCut; ++i) {
      assert(!newRow->isArcCut[i]);
      assert(!newRow->isArcCutReversed[i]);
    }
#endif
}

static void rigidFindStarNodes(MATRECNetworkDecomposition * dec, MATRECNetworkRowAddition * newRow,
                               const reduced_member_id toCheck){
    //4 cases:
    //Only a single edge; both head/tail are okay => network
    //All are adjacent to a single node, but do not have it as head or tail => not network
    //All are adjacent to a single node, and have it as head or tail => network
    //Not all are adjacent to a single node => check articulation nodes
    assert(newRow->reducedMembers[toCheck].numCutArcs > 0);//calling this function otherwise is nonsensical

    cut_arc_id cutArcIdx = newRow->reducedMembers[toCheck].firstCutArc;
    spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
    spqr_node head = findArcHead(dec, cutArc);
    spqr_node tail = findArcTail(dec, cutArc);

    bool reverse = findArcSign(dec,cutArc).reversed != newRow->cutArcs[cutArcIdx].arcReversed;
    spqr_node cutArcsHead = reverse ? tail : head;
    spqr_node cutArcsTail = reverse ? head : tail;

    if(newRow->reducedMembers[toCheck].numCutArcs == 1){
        newRow->reducedMembers[toCheck].articulationArc = cutArc;
        newRow->reducedMembers[toCheck].splitNode = cutArcsHead;
        newRow->reducedMembers[toCheck].otherNode = cutArcsTail;
        newRow->reducedMembers[toCheck].otherIsSource = true;
        newRow->reducedMembers[toCheck].allHaveCommonNode = true;
        return;// Only a single cut arc
    }

    spqr_node intersectionNodes[2] = {head,tail};

    while (cutArcIsValid(newRow->cutArcs[cutArcIdx].nextMember)) {
        cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
        cutArc = newRow->cutArcs[cutArcIdx].arc;
        head = findArcHead(dec, cutArc);
        tail = findArcTail(dec, cutArc);
        reverse = findArcSign(dec,cutArc).reversed != newRow->cutArcs[cutArcIdx].arcReversed;
        spqr_node effectiveHead = reverse ? tail : head;
        spqr_node effectiveTail = reverse ? head : tail;
        if(effectiveHead != cutArcsHead){
            cutArcsHead = SPQR_INVALID_NODE;
        }
        if(effectiveTail != cutArcsTail){
            cutArcsTail = SPQR_INVALID_NODE;
        }

        //intersection between intersectionNodes and head and tail
        for (int i = 0; i < 2; ++i) {
            if(intersectionNodes[i] != head && intersectionNodes[i] != tail){
                intersectionNodes[i] = SPQR_INVALID_NODE;
            }
        }
        if(SPQRnodeIsInvalid(intersectionNodes[0]) && SPQRnodeIsInvalid(intersectionNodes[1])){
            newRow->reducedMembers[toCheck].splitNode = SPQR_INVALID_NODE;
            newRow->reducedMembers[toCheck].allHaveCommonNode = false;
            return; //not all arcs are adjacent to a single node, need to check articulation nodes
        }
    }
    if(SPQRnodeIsInvalid(cutArcsHead) && SPQRnodeIsInvalid(cutArcsTail)){
        //All arcs adjacent to a single node, but not in same direction; not network
        newRow->remainsNetwork = false;
        newRow->reducedMembers[toCheck].type = TYPE_NOT_NETWORK;
        return;
    }
    bool headSplittable = SPQRnodeIsValid(cutArcsHead);
    //Check if the n arcs are in a n+1 degree node; if so, the other endpoint of this non split arc is also splittable
    //By virtue of the spanning tree, this arc must be a tree arc.
    spqr_node splitNode = headSplittable ? cutArcsHead : cutArcsTail;
    newRow->reducedMembers[toCheck].splitNode = splitNode;
    newRow->reducedMembers[toCheck].otherIsSource = headSplittable;
    newRow->reducedMembers[toCheck].allHaveCommonNode = true;
    if(newRow->reducedMembers[toCheck].numCutArcs == nodeDegree(dec,splitNode) - 1){
        spqr_arc firstNodeArc = getFirstNodeArc(dec,splitNode);
        spqr_arc neighbourArc = firstNodeArc;
        do{
            if(arcIsTree(dec,neighbourArc)){
                break;
            }
            neighbourArc = getNextNodeArc(dec,neighbourArc,splitNode);
        }while(neighbourArc != firstNodeArc);

        newRow->reducedMembers[toCheck].articulationArc = neighbourArc;
        spqr_arc arcHead = findArcHead(dec,neighbourArc);
        newRow->reducedMembers[toCheck].otherNode = arcHead == splitNode ? findArcTail(dec,neighbourArc) : arcHead;

    }
}

//TODO: remove MATREC_ERROR from below functions (until propagation function, basically) and refactor memory allocation
static MATREC_ERROR zeroOutColorsExceptNeighbourhood(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                                   const spqr_node articulationNode, const spqr_node startRemoveNode){
    COLOR_STATUS * neighbourColors;
    int degree = nodeDegree(dec,articulationNode);
    MATREC_CALL(MATRECallocBlockArray(dec->env,&neighbourColors,(size_t) degree));

    {
        int i = 0;
        spqr_arc artFirstArc = getFirstNodeArc(dec, articulationNode);
        spqr_arc artItArc = artFirstArc;
        do{
            spqr_node head = findArcHead(dec, artItArc);
            spqr_node tail = findArcTail(dec, artItArc);
            spqr_node otherNode = articulationNode == head ? tail : head;
            neighbourColors[i] = newRow->nodeColors[otherNode];
            i++;
            assert(i <= degree);
            artItArc = getNextNodeArc(dec,artItArc,articulationNode);
        }while(artItArc != artFirstArc);
    }
    zeroOutColors(dec,newRow,startRemoveNode);

    {
        int i = 0;
        spqr_arc artFirstArc = getFirstNodeArc(dec, articulationNode);
        spqr_arc artItArc = artFirstArc;
        do{
            spqr_node head = findArcHead(dec, artItArc);
            spqr_node tail = findArcTail(dec, artItArc);
            spqr_node otherNode = articulationNode == head ? tail : head;
            newRow->nodeColors[otherNode] = neighbourColors[i];
            i++;
            assert(i <= degree);
            artItArc = getNextNodeArc(dec,artItArc,articulationNode);
        }while(artItArc != artFirstArc);
    }

    MATRECfreeBlockArray(dec->env,&neighbourColors);
    return MATREC_OKAY;
}

static void intersectionOfAllPaths(MATRECNetworkDecomposition * dec, MATRECNetworkRowAddition *newRow,
                                   const reduced_member_id toCheck, int * const nodeNumPaths){
    int * intersectionPathDepth = newRow->intersectionPathDepth;
    spqr_node * intersectionPathParent = newRow->intersectionPathParent;

    //First do a dfs over the tree, storing all the tree-parents and depths for each node
    //TODO: maybe cache this tree and also update it so we can prevent this DFS call?


    //pick an arbitrary node as root; we just use the first cutArc here
    {
        spqr_node root = findArcHead(dec, newRow->cutArcs[newRow->reducedMembers[toCheck].firstCutArc].arc);
        DFSCallData *pathSearchCallStack = newRow->intersectionDFSData;

        assert(intersectionPathDepth[root] == -1);
        assert(intersectionPathParent[root] == SPQR_INVALID_NODE);

        int pathSearchCallStackSize = 0;

        intersectionPathDepth[root] = 0;
        intersectionPathParent[root] = SPQR_INVALID_NODE;

        pathSearchCallStack[0].node = root;
        pathSearchCallStack[0].nodeArc = getFirstNodeArc(dec, root);
        pathSearchCallStackSize++;
        while (pathSearchCallStackSize > 0) {
            assert(pathSearchCallStackSize <= newRow->memIntersectionDFSData);
            DFSCallData *dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
            //cannot be a tree arc which is its parent
            if (arcIsTree(dec, dfsData->nodeArc) &&
                (pathSearchCallStackSize <= 1 ||
                 dfsData->nodeArc != pathSearchCallStack[pathSearchCallStackSize - 2].nodeArc)) {
                spqr_node head = findArcHeadNoCompression(dec, dfsData->nodeArc);
                spqr_node tail = findArcTailNoCompression(dec, dfsData->nodeArc);
                spqr_node other = head == dfsData->node ? tail : head;
                assert(other != dfsData->node);

                //We go up a level: add new node to the call stack
                pathSearchCallStack[pathSearchCallStackSize].node = other;
                pathSearchCallStack[pathSearchCallStackSize].nodeArc = getFirstNodeArc(dec, other);
                //Every time a new node is discovered/added, we update its parent and depth information
                assert(intersectionPathDepth[other] == -1);
                assert(intersectionPathParent[other] == SPQR_INVALID_NODE);
                intersectionPathParent[other] = dfsData->node;
                intersectionPathDepth[other] = pathSearchCallStackSize;
                ++pathSearchCallStackSize;
                continue;
            }
            do {
                dfsData->nodeArc = getNextNodeArc(dec, dfsData->nodeArc, dfsData->node);
                if (dfsData->nodeArc == getFirstNodeArc(dec, dfsData->node)) {
                    --pathSearchCallStackSize;
                    dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
                } else {
                    break;
                }
            } while (pathSearchCallStackSize > 0);
        }
    }

    //For each cut arc, trace back both ends until they meet
    cut_arc_id cutArc = newRow->reducedMembers[toCheck].firstCutArc;
    do{
        spqr_arc arc = newRow->cutArcs[cutArc].arc;
        cutArc = newRow->cutArcs[cutArc].nextMember;

        //Iteratively jump up to the parents until they reach a common parent
        spqr_node source = findArcHead(dec, arc);
        spqr_node target = findArcTail(dec, arc);
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

    }while (cutArcIsValid(cutArc));

}

static void addArticulationNode(MATRECNetworkRowAddition *newRow, spqr_node articulationNode){
#ifndef NDEBUG
    for (int i = 0; i < newRow->numArticulationNodes; ++i) {
        assert(newRow->articulationNodes[i] != articulationNode);
    }
#endif
    newRow->articulationNodes[newRow->numArticulationNodes] = articulationNode;
    ++newRow->numArticulationNodes;
}
static void articulationPoints(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition * newRow, ArticulationNodeInformation *nodeInfo, reduced_member_id reducedMember){
    const bool * arcRemoved = newRow->isArcCut;

    int rootChildren = 0;
    spqr_node root_node = findArcHead(dec, getFirstMemberArc(dec, newRow->reducedMembers[reducedMember].member));;

    ArticulationPointCallStack * callStack = newRow->artDFSData;

    int depth = 0;
    int time = 1;

    callStack[depth].arc = getFirstNodeArc(dec,root_node);
    callStack[depth].node = root_node;
    callStack[depth].parent = SPQR_INVALID_NODE;
    callStack[depth].isAP = false;

    nodeInfo[root_node].low = time;
    nodeInfo[root_node].discoveryTime = time;

    while(depth >= 0){
        if(!arcRemoved[callStack[depth].arc]){
            spqr_node node = callStack[depth].node;
            spqr_node head = findArcHead(dec, callStack[depth].arc);
            spqr_node tail = findArcTail(dec, callStack[depth].arc);
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
                    callStack[depth].arc = getFirstNodeArc(dec,otherNode);
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
            callStack[depth].arc = getNextNodeArc(dec,callStack[depth].arc,callStack[depth].node);
            if(callStack[depth].arc != getFirstNodeArc(dec,callStack[depth].node)) break;
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

static void rigidConnectedColoringRecursive(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition * newRow, spqr_node articulationNode,
                                            spqr_node firstProcessNode, COLOR_STATUS firstColor,
                                            bool *isGood){
    const bool * isArcCut = newRow->isArcCut;
    COLOR_STATUS * nodeColors = newRow->nodeColors;
    ColorDFSCallData * data = newRow->colorDFSData;

    data[0].node = firstProcessNode;
    data[0].arc = getFirstNodeArc(dec,firstProcessNode);
    newRow->nodeColors[firstProcessNode] = firstColor;

    int depth = 0;
    while(depth >= 0){
        assert(depth < newRow->memColorDFSData);
        assert(newRow->nodeColors[articulationNode] == UNCOLORED);

        ColorDFSCallData * callData = &data[depth];
        spqr_node head = findArcHead(dec, callData->arc);
        spqr_node tail = findArcTail(dec, callData->arc);
        spqr_node otherNode = callData->node == head ? tail : head;
        COLOR_STATUS currentColor = nodeColors[callData->node];
        COLOR_STATUS otherColor = nodeColors[otherNode];
        //Checks the direction of the arc; in the rest of the algorithm, we just need to check partition
        if(isArcCut[callData->arc] && currentColor != otherColor){
            bool otherIsTail = callData->node == head;
            bool arcReversed = findArcSign(dec,callData->arc).reversed != newRow->isArcCutReversed[callData->arc];
            bool good = (currentColor == COLOR_SOURCE) == (otherIsTail == arcReversed);
            if(!good){
                *isGood = false;
                break;
            }
        }
        if(otherNode != articulationNode){
            if(otherColor == UNCOLORED){
                if(isArcCut[callData->arc]){
                    nodeColors[otherNode] = currentColor == COLOR_SOURCE ? COLOR_SINK : COLOR_SOURCE; //reverse the colors
                }else{
                    nodeColors[otherNode] = currentColor;
                }
                callData->arc = getNextNodeArc(dec,callData->arc,callData->node);

                depth++;
                assert(depth < newRow->memColorDFSData);
                data[depth].node = otherNode;
                data[depth].arc = getFirstNodeArc(dec,otherNode);
                continue;
            }
            if(isArcCut[callData->arc] != (currentColor != otherColor)){
                *isGood = false;
                break;
            }
        }
        callData->arc = getNextNodeArc(dec,callData->arc,callData->node);
        while(depth >= 0 && data[depth].arc == getFirstNodeArc(dec,data[depth].node)){
            --depth;
        }
    }
}

static void rigidConnectedColoring(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                   const reduced_member_id reducedMember, const spqr_node node, bool * const isGood){

    //we should only perform this function if there's more than one cut arc
    assert(newRow->reducedMembers[reducedMember].numCutArcs > 1);
#ifndef NDEBUG
    {

        spqr_member member = newRow->reducedMembers[reducedMember].member;
        spqr_arc firstArc = getFirstMemberArc(dec, member);
        spqr_arc memberArc = firstArc;
        do{
            assert(newRow->nodeColors[findArcHeadNoCompression(dec,memberArc)] == UNCOLORED);
            assert(newRow->nodeColors[findArcTailNoCompression(dec,memberArc)] == UNCOLORED);
            memberArc = getNextMemberArc(dec,memberArc);
        }while(firstArc != memberArc);
    }
#endif

    spqr_node firstProcessNode;
    COLOR_STATUS firstColor;
    {
        cut_arc_id cutArc = newRow->reducedMembers[reducedMember].firstCutArc;
        spqr_arc arc = newRow->cutArcs[cutArc].arc;
        assert(SPQRarcIsValid(arc));
        spqr_node head = findArcHead(dec, arc);
        spqr_node tail = findArcTail(dec, arc);
        if(findArcSign(dec,arc).reversed != newRow->cutArcs[cutArc].arcReversed){
            spqr_node temp = head;
            head = tail;
            tail = temp;
        }
        if(tail != node){
            firstProcessNode = tail;
            firstColor = COLOR_SOURCE;
        }else{
            assert(head != node);
            firstProcessNode = head;
            firstColor = COLOR_SINK;
        }
    }
    assert(SPQRnodeIsValid(firstProcessNode) && firstProcessNode != node);
    *isGood = true;
    newRow->reducedMembers[reducedMember].coloredNode = firstProcessNode;
    rigidConnectedColoringRecursive(dec,newRow,node,firstProcessNode,firstColor,isGood);

    // Need to zero all colors for next attempts if we failed
    if(!(*isGood)){
        zeroOutColors(dec,newRow,firstProcessNode);
        newRow->reducedMembers[reducedMember].coloredNode = SPQR_INVALID_NODE;
    }else{
        //Otherwise, we zero out all colors but the ones which we need
        zeroOutColorsExceptNeighbourhood(dec,newRow,node, firstProcessNode);
        newRow->reducedMembers[reducedMember].coloredNode = node;
    }
}

static spqr_node checkNeighbourColoringArticulationNode(MATRECNetworkDecomposition * dec, MATRECNetworkRowAddition *newRow,
                                                        const spqr_node articulationNode, spqr_arc * const adjacentSplittingArc){
    spqr_node firstSideCandidate = SPQR_INVALID_NODE;
    spqr_node secondSideCandidate = SPQR_INVALID_NODE;
    spqr_arc firstSideArc = SPQR_INVALID_ARC;
    spqr_arc secondSideArc = SPQR_INVALID_ARC;
    int numFirstSide = 0;
    int numSecondSide = 0;

    spqr_arc firstArc = getFirstNodeArc(dec, articulationNode);
    spqr_arc moveArc = firstArc;
    do{
        spqr_node head = findArcHead(dec, moveArc);
        spqr_node tail = findArcTail(dec, moveArc);
        spqr_node otherNode = articulationNode == head ? tail : head;
        assert(newRow->nodeColors[otherNode] != UNCOLORED);
        if((newRow->nodeColors[otherNode] == COLOR_SOURCE) != newRow->isArcCut[moveArc] ){
            if(numFirstSide == 0 && arcIsTree(dec,moveArc)){
                firstSideCandidate = otherNode;
                firstSideArc = moveArc;
            }else if (numFirstSide > 0){
                firstSideCandidate = SPQR_INVALID_NODE;
            }
            ++numFirstSide;
        }else{
            if(numSecondSide == 0 && arcIsTree(dec,moveArc)){
                secondSideCandidate = otherNode;
                secondSideArc = moveArc;
            }else if (numSecondSide > 0){
                secondSideCandidate = SPQR_INVALID_NODE;
            }
            ++numSecondSide;
        }
        moveArc = getNextNodeArc(dec,moveArc,articulationNode);
    }while(moveArc != firstArc);

    if(numFirstSide == 1){
        *adjacentSplittingArc = firstSideArc;
        return firstSideCandidate;
    }else if (numSecondSide == 1){
        *adjacentSplittingArc = secondSideArc;
        return secondSideCandidate;
    }
    return SPQR_INVALID_NODE;
}


static void rigidGetSplittableArticulationPointsOnPath(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                                       const reduced_member_id toCheck,
                                                       spqr_node firstNode, spqr_node secondNode){
    assert(newRow->reducedMembers[toCheck].numCutArcs > 1);

    int totalNumNodes = getNumNodes(dec);
    int * nodeNumPaths = newRow->crossingPathCount;

    for (int i = 0; i < totalNumNodes; ++i) {
        nodeNumPaths[i] = 0;
    }

    intersectionOfAllPaths(dec,newRow,toCheck,nodeNumPaths);

    newRow->numArticulationNodes = 0;

    ArticulationNodeInformation * artNodeInfo = newRow->articulationNodeSearchInfo;
    //TODO: ugly; we clean up over all decomposition nodes for every component
    //clean up can not easily be done in the search, unfortunately
    for (int i = 0; i < totalNumNodes; ++i) {
        artNodeInfo[i].low = 0 ;
        artNodeInfo[i].discoveryTime = 0;
    }

    articulationPoints(dec,newRow,artNodeInfo,toCheck);

    int numCutArcs = newRow->reducedMembers[toCheck].numCutArcs;
    for (int i = 0; i < newRow->numArticulationNodes; i++) {
        spqr_node articulationNode = newRow->articulationNodes[i];
        assert(nodeIsRepresentative(dec, articulationNode));
        bool isOnPath = nodeNumPaths[articulationNode] == numCutArcs;
        if (isOnPath && (
                (SPQRnodeIsInvalid(firstNode) &&SPQRnodeIsInvalid(secondNode))
                || articulationNode == firstNode || articulationNode == secondNode )){
            bool isGood = true;
            rigidConnectedColoring(dec, newRow, toCheck, articulationNode, &isGood);
            if(!isGood){
                continue;
            }
            newRow->reducedMembers[toCheck].splitNode = articulationNode;

            //Once we have found one node, we can (possibly) find another by looking at the coloring of the neighbourhood of it.

            spqr_arc adjacentSplittingArc = SPQR_INVALID_ARC;
            spqr_node adjacentSplittingNode = checkNeighbourColoringArticulationNode(dec, newRow, articulationNode,
                                                                                     &adjacentSplittingArc);

            //If the neighbour-coloring works, we still need to check that the adjacent node
            //is also an articulation node
            if(SPQRnodeIsValid(adjacentSplittingNode) && ((SPQRnodeIsInvalid(firstNode) && SPQRnodeIsInvalid(secondNode))
               || adjacentSplittingNode == firstNode || adjacentSplittingNode == secondNode)){
                bool isArticulationNode = false;
                for (int j = 0; j < newRow->numArticulationNodes; ++j) {
                    if (newRow->articulationNodes[j] == adjacentSplittingNode) {
                        isArticulationNode = true;
                        break;
                    }
                }
                if (isArticulationNode) {
                    newRow->reducedMembers[toCheck].articulationArc = adjacentSplittingArc;
                    newRow->reducedMembers[toCheck].otherNode = adjacentSplittingNode;
                    newRow->reducedMembers[toCheck].otherIsSource = newRow->nodeColors[adjacentSplittingNode] == COLOR_SOURCE;

                    //Cleaning up the colors
                    {
                        spqr_arc firstNodeArc = getFirstNodeArc(dec, articulationNode);
                        spqr_arc itArc = firstNodeArc;
                        do {
                            spqr_node head = findArcHead(dec, itArc);
                            spqr_node tail = findArcTail(dec, itArc);
                            spqr_node otherNode = articulationNode == head ? tail : head;
                            newRow->nodeColors[otherNode] = UNCOLORED;
                            itArc = getNextNodeArc(dec, itArc, articulationNode);
                        } while (itArc != firstNodeArc);

                    }
                }
            }
            return;
        }
    }

    newRow->reducedMembers[toCheck].type = TYPE_NOT_NETWORK;
    newRow->remainsNetwork = false;
}

static void determineParallelType(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                  const reduced_member_id toCheckMember, const spqr_arc markerToOther,
                                  const reduced_member_id otherMember, const spqr_arc markerToCheck){

    MATRECRowReducedMember * member = &newRow->reducedMembers[toCheckMember];
    assert(cutArcIsValid(member->firstCutArc));

    bool good = true;
    bool isReversed = true;
    int countedCutArcs = 0;
    for (cut_arc_id cutArc = member->firstCutArc; cutArcIsValid(cutArc);
         cutArc = newRow->cutArcs[cutArc].nextMember){
        spqr_arc arc = newRow->cutArcs[cutArc].arc;
        bool arcIsReversed = arcIsReversedNonRigid(dec,arc) != newRow->cutArcs[cutArc].arcReversed;
        if(countedCutArcs == 0){
            isReversed = arcIsReversed;
        }else if(arcIsReversed != isReversed){
            good = false;
            break;
        }
        ++countedCutArcs;
    }
    if(!good){
        member->type = TYPE_NOT_NETWORK;
        newRow->remainsNetwork = false;
        return;
    }
    //we can only propagate if the marker arc is a tree arc and all other arcs are cut
    if(!arcIsTree(dec,markerToOther) ||
       countedCutArcs != (getNumMemberArcs(dec,newRow->reducedMembers[toCheckMember].member) - 1)){
        //In all other cases, the bond can be split so that the result will be okay!
        newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
    }else{
        bool markerIsReversed = arcIsReversedNonRigid(dec,markerToOther);
        createCutArc(dec,newRow,markerToCheck,otherMember, markerIsReversed != isReversed);
        newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
    }

}

static void determineSeriesType(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                const reduced_member_id toCheckMember, const spqr_arc markerToOther,
                                const reduced_member_id otherMember, const spqr_arc markerToCheck){
    //Propagation only calls this function if the arc is tree already, so we do not check it here.
    assert(arcIsTree(dec,markerToOther));
    assert(newRow->reducedMembers[toCheckMember].numCutArcs == 1);
    assert(cutArcIsValid(newRow->reducedMembers[toCheckMember].firstCutArc));
    spqr_arc cutArc = newRow->reducedMembers[toCheckMember].firstCutArc;
    spqr_arc arc = newRow->cutArcs[cutArc].arc;
    newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
    createCutArc(dec,newRow,markerToCheck,otherMember,
                 (arcIsReversedNonRigid(dec,arc) == arcIsReversedNonRigid(dec,markerToOther)) != newRow->cutArcs[cutArc].arcReversed);
}

static void determineRigidType(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                               const reduced_member_id toCheckMember, const spqr_arc markerToOther,
                               const reduced_member_id otherMember, const spqr_arc markerToCheck){
    assert(newRow->reducedMembers[toCheckMember].numCutArcs > 0);//Checking for propagation only makes sense if there is at least one cut arc

    rigidFindStarNodes(dec,newRow,toCheckMember);
    if(!newRow->remainsNetwork){
        return;
    }
    spqr_node markerHead = findArcHead(dec,markerToOther);
    spqr_node markerTail = findArcTail(dec,markerToOther);
    if(findArcSign(dec,markerToOther).reversed){
        spqr_node temp = markerHead;
        markerHead = markerTail;
        markerTail = temp;
    }
    if(SPQRnodeIsInvalid(newRow->reducedMembers[toCheckMember].splitNode)){
       //not a star => attempt to find splittable nodes using articulation node algorithms
       //We save some work by telling the methods that only the marker nodes should be checked
       rigidGetSplittableArticulationPointsOnPath(dec,newRow,toCheckMember,markerHead,markerTail);
    }
    if(!newRow->remainsNetwork){
        return;
    }


    if(SPQRarcIsValid(newRow->reducedMembers[toCheckMember].articulationArc) &&
            newRow->reducedMembers[toCheckMember].articulationArc == markerToOther){
        newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
        bool reverse =(markerHead == newRow->reducedMembers[toCheckMember].splitNode) !=
                newRow->reducedMembers[toCheckMember].otherIsSource;

        createCutArc(dec,newRow,markerToCheck,otherMember,reverse);
    }else if(newRow->reducedMembers[toCheckMember].splitNode == markerHead ||
            newRow->reducedMembers[toCheckMember].splitNode == markerTail){
        newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
    }else if(SPQRarcIsValid(newRow->reducedMembers[toCheckMember].articulationArc) &&
            (newRow->reducedMembers[toCheckMember].otherNode == markerHead ||
            newRow->reducedMembers[toCheckMember].otherNode == markerTail)
    ){
        newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
    }else{
        //Found source or sinks, but not adjacent to the marker
        newRow->reducedMembers[toCheckMember].type = TYPE_NOT_NETWORK;
        newRow->remainsNetwork = false;
    }

}
static RowReducedMemberType determineType(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                          const reduced_member_id toCheckMember, const spqr_arc markerToOther,
                                          const reduced_member_id otherMember, const spqr_arc markerToCheck){
    assert(newRow->reducedMembers[toCheckMember].type == TYPE_UNDETERMINED);
    switch (getMemberType(dec,newRow->reducedMembers[toCheckMember].member)) {
        case SPQR_MEMBERTYPE_RIGID:
        {
            determineRigidType(dec,newRow,toCheckMember,markerToOther,otherMember,markerToCheck);
            break;

        }
        case SPQR_MEMBERTYPE_PARALLEL:
        {
            determineParallelType(dec,newRow,toCheckMember,markerToOther,otherMember,markerToCheck);
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:
        {
            determineSeriesType(dec,newRow,toCheckMember,markerToOther,otherMember,markerToCheck);
            break;
        }
        default:
            assert(false);
            newRow->reducedMembers[toCheckMember].type = TYPE_NOT_NETWORK;
    }

    return newRow->reducedMembers[toCheckMember].type;
}

static void propagateComponents(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow){
    int leafArrayIndex = 0;

    reduced_member_id leaf;
    reduced_member_id next;

    while(leafArrayIndex != newRow->numLeafMembers){
        leaf = newRow->leafMembers[leafArrayIndex];
        //next is invalid if the member is not in the reduced decomposition.
        next = newRow->reducedMembers[leaf].parent;
        spqr_arc parentMarker = markerToParent(dec, newRow->reducedMembers[leaf].member);
        if(next != INVALID_REDUCED_MEMBER && arcIsTree(dec,parentMarker)){
            assert(reducedMemberIsValid(next));
            assert(SPQRarcIsValid(parentMarker));
            RowReducedMemberType type = determineType(dec,newRow,leaf,parentMarker,next,markerOfParent(dec,newRow->reducedMembers[leaf].member));
            if(type == TYPE_PROPAGATED){
                ++newRow->reducedMembers[next].numPropagatedChildren;
                if(newRow->reducedMembers[next].numPropagatedChildren == newRow->reducedMembers[next].numChildren){
                    newRow->leafMembers[leafArrayIndex] = next;
                }else{
                    ++leafArrayIndex;
                }
            }else if(type == TYPE_NOT_NETWORK){
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
                spqr_arc markerToChild = SPQR_INVALID_ARC;
                for (children_idx i = newRow->reducedMembers[root].firstChild; i < newRow->reducedMembers[root].firstChild + newRow->reducedMembers[root].numChildren; ++i) {
                    reduced_member_id childReduced = newRow->childrenStorage[i];
                    if(newRow->reducedMembers[childReduced].type != TYPE_PROPAGATED){
                        child = childReduced;
                        markerToChild = markerOfParent(dec,newRow->reducedMembers[child].member);
                        break;
                    }
                }
                assert(SPQRmemberIsValid(newRow->reducedMembers[child].member));
                assert(SPQRarcIsValid(markerToChild));
                if(!arcIsTree(dec,markerToChild)){
                    break;
                }
                RowReducedMemberType type = determineType(dec,newRow,root,markerToChild,child,markerToParent(dec,newRow->reducedMembers[child].member));
                if(type == TYPE_PROPAGATED){
                    root = child;
                }else if(type == TYPE_NOT_NETWORK){
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

typedef struct {
    spqr_node first;
    spqr_node second;
} Nodes;
static Nodes
rigidDetermineCandidateNodesFromAdjacentComponents(MATRECNetworkDecomposition *dec,
                                                   MATRECNetworkRowAddition *newRow, const reduced_member_id toCheck) {

    Nodes pair;
    pair.first = SPQR_INVALID_NODE;
    pair.second = SPQR_INVALID_NODE;

    //take union of children's arcs nodes to find one or two candidates
    for (children_idx i = newRow->reducedMembers[toCheck].firstChild;
         i < newRow->reducedMembers[toCheck].firstChild + newRow->reducedMembers[toCheck].numChildren; ++i) {
        reduced_member_id reducedChild = newRow->childrenStorage[i];
        if (newRow->reducedMembers[reducedChild].type != TYPE_PROPAGATED) {
            spqr_arc arc = markerOfParent(dec, newRow->reducedMembers[reducedChild].member);
            spqr_node head = findArcHead(dec, arc);
            spqr_node tail = findArcTail(dec, arc);
            if(SPQRnodeIsInvalid(pair.first) && SPQRnodeIsInvalid(pair.second)){
                pair.first = head;
                pair.second = tail;
            }else{
                if(pair.first != head && pair.first != tail){
                    pair.first = SPQR_INVALID_NODE;
                }
                if(pair.second != head && pair.second != tail){
                    pair.second = SPQR_INVALID_NODE;
                }
            }
            if (SPQRnodeIsInvalid(pair.first) && SPQRnodeIsInvalid(pair.second)) {
                return pair;
            }
        }
    }
    if (reducedMemberIsValid(newRow->reducedMembers[toCheck].parent) &&
        newRow->reducedMembers[newRow->reducedMembers[toCheck].parent].type != TYPE_PROPAGATED) {

        spqr_arc arc = markerToParent(dec, newRow->reducedMembers[toCheck].member);
        spqr_node head = findArcHead(dec, arc);
        spqr_node tail = findArcTail(dec, arc);
        if(SPQRnodeIsInvalid(pair.first) && SPQRnodeIsInvalid(pair.second)){
            pair.first = head;
            pair.second = tail;
        }else{
            if(pair.first != head && pair.first != tail){
                pair.first = SPQR_INVALID_NODE;
            }
            if(pair.second != head && pair.second != tail){
                pair.second = SPQR_INVALID_NODE;
            }
        }
    }
    return pair;
}

static MATREC_ERROR allocateTreeSearchMemory(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow){
    int necessarySpace = newRow->numReducedMembers;
    if( necessarySpace > newRow->memMergeTreeCallData ){
        newRow->memMergeTreeCallData = max(2*newRow->memMergeTreeCallData,necessarySpace);
        MATREC_CALL(MATRECreallocBlockArray(dec->env,&newRow->mergeTreeCallData,(size_t) newRow->memMergeTreeCallData));
    }
    return MATREC_OKAY;
}

static void determineSingleRowRigidType(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                        reduced_member_id reducedMember){
    assert(newRow->reducedMembers[reducedMember].numCutArcs > 0);//Checking for propagation only makes sense if there is at least one cut arc
    rigidFindStarNodes(dec,newRow,reducedMember);
    if(!newRow->remainsNetwork){
        return;
    }

    if(SPQRnodeIsInvalid(newRow->reducedMembers[reducedMember].splitNode)){
        //not a star => attempt to find splittable nodes using articulation node algorithms
        rigidGetSplittableArticulationPointsOnPath(dec,newRow,reducedMember,SPQR_INVALID_NODE,SPQR_INVALID_NODE);
    }
    if(SPQRnodeIsValid(newRow->reducedMembers[reducedMember].splitNode)){
        newRow->reducedMembers[reducedMember].type = TYPE_MERGED;
    }else{
        newRow->reducedMembers[reducedMember].type = TYPE_NOT_NETWORK;
        newRow->remainsNetwork = false;
    }

}
static void determineSingleParallelType(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                        reduced_member_id reducedMember){
    MATRECRowReducedMember * redMember = &newRow->reducedMembers[reducedMember];
    assert(cutArcIsValid(redMember->firstCutArc));

    bool good = true;
    bool isReversed = true;
    int countedCutArcs = 0;
    for (cut_arc_id cutArc = redMember->firstCutArc; cutArcIsValid(cutArc);
         cutArc = newRow->cutArcs[cutArc].nextMember){
        spqr_arc arc = newRow->cutArcs[cutArc].arc;
        bool arcIsReversed = arcIsReversedNonRigid(dec,arc) != newRow->cutArcs[cutArc].arcReversed;
        if(countedCutArcs == 0){
            isReversed = arcIsReversed;
        }else if(arcIsReversed != isReversed){
            good = false;
            break;
        }
        ++countedCutArcs;
    }
    if(!good){
        redMember->type = TYPE_NOT_NETWORK;
        newRow->remainsNetwork = false;
    }else{
        redMember->type = TYPE_MERGED;
    }
}
static void determineSingleSeriesType(MATRECNetworkRowAddition *newRow,
                                        reduced_member_id reducedMember){
    assert(newRow->reducedMembers[reducedMember].numCutArcs == 1);
    assert(cutArcIsValid(newRow->reducedMembers[reducedMember].firstCutArc));
    newRow->reducedMembers[reducedMember].type = TYPE_MERGED;
}

static spqr_node determineAndColorSplitNode(MATRECNetworkDecomposition * dec, MATRECNetworkRowAddition * newRow,
                                       reduced_member_id id, spqr_node candidateOne, spqr_node candidateTwo){
    if(SPQRnodeIsInvalid(newRow->reducedMembers[id].splitNode)){
        return SPQR_INVALID_NODE;
    }
    spqr_node splitNode = newRow->reducedMembers[id].splitNode;
    if(splitNode == candidateOne || splitNode == candidateTwo){
        if(newRow->reducedMembers[id].allHaveCommonNode){
            spqr_arc firstNodeArc = getFirstNodeArc(dec,splitNode);
            spqr_arc iterArc = firstNodeArc;
            COLOR_STATUS color = newRow->reducedMembers[id].otherIsSource ? COLOR_SOURCE : COLOR_SINK;
            do{
                spqr_node head = findArcHead(dec,iterArc);
                spqr_node other = head == splitNode ? findArcTail(dec,iterArc) : head;
                newRow->nodeColors[other] = color;
                iterArc = getNextNodeArc(dec,iterArc,splitNode);
            }while(iterArc != firstNodeArc);
            newRow->reducedMembers[id].coloredNode = splitNode;
        }
        return splitNode;
    }
    splitNode = newRow->reducedMembers[id].otherNode;
    if(SPQRnodeIsInvalid(splitNode) || (splitNode != candidateOne && splitNode != candidateTwo)){
        return SPQR_INVALID_NODE;
    }
    assert(SPQRarcIsValid(newRow->reducedMembers[id].articulationArc));
    if(newRow->reducedMembers[id].allHaveCommonNode){
        spqr_arc firstNodeArc = getFirstNodeArc(dec,splitNode);
        spqr_arc iterArc = firstNodeArc;
        COLOR_STATUS color;
        if(newRow->reducedMembers[id].numCutArcs == 1){
            color = newRow->reducedMembers[id].otherIsSource ? COLOR_SINK : COLOR_SOURCE;
        }else{
            color = newRow->reducedMembers[id].otherIsSource ? COLOR_SOURCE : COLOR_SINK;
        }
        do{
            spqr_node head = findArcHead(dec,iterArc);
            spqr_node other = head == splitNode ? findArcTail(dec,iterArc) : head;
            newRow->nodeColors[other] = color;
            iterArc = getNextNodeArc(dec,iterArc,splitNode);
        }while(iterArc != firstNodeArc);
        newRow->nodeColors[newRow->reducedMembers[id].splitNode] = newRow->reducedMembers[id].otherIsSource ? COLOR_SINK : COLOR_SOURCE;

    }else{
        COLOR_STATUS splitColor = newRow->nodeColors[splitNode];

        spqr_arc firstNodeArc = getFirstNodeArc(dec,splitNode);
        spqr_arc iterArc = firstNodeArc;
        do{
            spqr_node head = findArcHead(dec,iterArc);
            spqr_node other = head == splitNode ? findArcTail(dec,iterArc) : head;
            newRow->nodeColors[other] = splitColor;
            iterArc = getNextNodeArc(dec,iterArc,splitNode);
        }while(iterArc != firstNodeArc);
        newRow->nodeColors[newRow->reducedMembers[id].splitNode] = splitColor == COLOR_SOURCE ? COLOR_SINK : COLOR_SOURCE;
        newRow->nodeColors[splitNode] = UNCOLORED;
    }

    newRow->reducedMembers[id].coloredNode = splitNode;
    return splitNode;
}
static void determineSplitTypeFirstLeaf(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                        reduced_member_id reducedId){
    spqr_member member = newRow->reducedMembers[reducedId].member;
    SPQRMemberType type = getMemberType(dec,member);
    assert(type == SPQR_MEMBERTYPE_PARALLEL || type == SPQR_MEMBERTYPE_RIGID);
    assert(cutArcIsValid(newRow->reducedMembers[reducedId].firstCutArc));
    MATRECRowReducedMember * redMember = &newRow->reducedMembers[reducedId];

    if(type == SPQR_MEMBERTYPE_PARALLEL){
        //TODO: duplicate-ish

        bool good = true;
        bool isReversed = true;
        int countedCutArcs = 0;
        for (cut_arc_id cutArc = redMember->firstCutArc; cutArcIsValid(cutArc);
             cutArc = newRow->cutArcs[cutArc].nextMember){
            spqr_arc arc = newRow->cutArcs[cutArc].arc;
            bool arcIsReversed = arcIsReversedNonRigid(dec,arc) != newRow->cutArcs[cutArc].arcReversed;
            if(countedCutArcs == 0){
                isReversed = arcIsReversed;
            }else if(arcIsReversed != isReversed){
                good = false;
                break;
            }
            ++countedCutArcs;
        }
        if(!good){
            redMember->type = TYPE_NOT_NETWORK;
            newRow->remainsNetwork = false;
        }else{
            spqr_arc marker = markerToParent(dec,member);
            redMember->type = TYPE_MERGED;
            redMember->splitArc = marker;
            redMember->splitHead = true;
            redMember->otherIsSource = arcIsReversedNonRigid(dec,marker) == isReversed;
        }
        return;
    }
    assert(type == SPQR_MEMBERTYPE_RIGID);

    spqr_arc marker = markerToParent(dec,member);
    spqr_node markerHead = findArcHead(dec,marker);
    spqr_node markerTail = findArcTail(dec,marker);
    if(findArcSign(dec,marker).reversed){
        spqr_node temp = markerHead;
        markerHead = markerTail;
        markerTail = temp;
    }

    if(!SPQRnodeIsValid(newRow->reducedMembers[reducedId].splitNode)){
        assert(newRow->reducedMembers[reducedId].numCutArcs > 0);//Checking for propagation only makes sense if there is at least one cut arc

        rigidFindStarNodes(dec,newRow,reducedId);
        if(!newRow->remainsNetwork){
            return;
        }

        if(SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode)){
            //not a star => attempt to find splittable nodes using articulation node algorithms
            //We save some work by telling the methods that only the marker nodes should be checked
            rigidGetSplittableArticulationPointsOnPath(dec,newRow,reducedId,markerHead,markerTail);
        }
        if(!newRow->remainsNetwork){
            return;
        }
        if(SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode)){
            redMember->type = TYPE_NOT_NETWORK;
            newRow->remainsNetwork = false;
            return;
        }

    }

    spqr_node splitNode = newRow->reducedMembers[reducedId].splitNode;
    spqr_node otherNode = newRow->reducedMembers[reducedId].otherNode;
    //We cannot have both splittable (should have been propagated)
    assert(!((splitNode == markerTail || splitNode == markerHead)
             && (otherNode == markerTail || otherNode == markerHead)));

    splitNode = determineAndColorSplitNode(dec,newRow,reducedId,markerHead,markerTail);
    if(SPQRnodeIsInvalid(splitNode)){
        redMember->type = TYPE_NOT_NETWORK;
        newRow->remainsNetwork = false;
        return;
    }
    assert(splitNode == markerHead || splitNode == markerTail);

    newRow->reducedMembers[reducedId].splitNode = splitNode;
    newRow->reducedMembers[reducedId].willBeReversed = false;
    redMember->type = TYPE_MERGED;

}
typedef struct{
    bool headSplit;
    bool otherIsSource;
} SplitOrientation;

static SplitOrientation getRelativeOrientationRigid(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                                       reduced_member_id reducedId, spqr_arc arcToNext){
    assert(findArcMemberNoCompression(dec,arcToNext) == newRow->reducedMembers[reducedId].member);
    assert(SPQRnodeIsValid(newRow->reducedMembers[reducedId].splitNode));

    //TODO: can probably simplify, but we only do this after we are sure all the cases are correct
    SplitOrientation orientation;
    if(newRow->reducedMembers[reducedId].numCutArcs == 0){
        spqr_node splitNode = newRow->reducedMembers[reducedId].splitNode;
        spqr_node head = findEffectiveArcHead(dec,arcToNext);

        assert(head == splitNode || splitNode == findEffectiveArcTailNoCompression(dec,arcToNext));

        orientation.headSplit = newRow->reducedMembers[reducedId].willBeReversed == (head != splitNode);
        orientation.otherIsSource = newRow->reducedMembers[reducedId].otherIsSource;
        return orientation;
    }
    spqr_node splitNode = newRow->reducedMembers[reducedId].splitNode;
    spqr_node arcHead = findArcHead(dec,arcToNext);
    spqr_node arcTail = findArcTail(dec,arcToNext);
    if(findArcSign(dec,arcToNext).reversed){
        spqr_node temp = arcHead;
        arcHead = arcTail;
        arcTail = temp;
    }
    assert(arcHead == splitNode || arcTail == splitNode);
    spqr_node other = arcHead == splitNode ? arcTail : arcHead;

    assert(newRow->nodeColors[other] == COLOR_SOURCE || newRow->nodeColors[other] == COLOR_SINK);

    if(newRow->reducedMembers[reducedId].willBeReversed){
        orientation.headSplit = arcHead != splitNode;
        orientation.otherIsSource = newRow->nodeColors[other] != COLOR_SOURCE;
    }else{
        orientation.headSplit = arcHead == splitNode;
        orientation.otherIsSource =  newRow->nodeColors[other] == COLOR_SOURCE;
    }
    return orientation;
}

static SplitOrientation getRelativeOrientationParallel(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                               reduced_member_id reducedId, spqr_arc arcToNext){
    assert(findArcMemberNoCompression(dec,arcToNext) == newRow->reducedMembers[reducedId].member);
    assert(SPQRarcIsValid(newRow->reducedMembers[reducedId].splitArc) && SPQRarcIsValid(arcToNext));
    SplitOrientation orientation;
    orientation.otherIsSource = newRow->reducedMembers[reducedId].otherIsSource;
    if(arcIsReversedNonRigid(dec,arcToNext) == arcIsReversedNonRigid(dec,newRow->reducedMembers[reducedId].splitArc)){
        orientation.headSplit = newRow->reducedMembers[reducedId].splitHead;
    }else{
        orientation.headSplit = !newRow->reducedMembers[reducedId].splitHead;
    }
    return orientation;
}
static SplitOrientation getRelativeOrientationSeries(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                                       reduced_member_id reducedId, spqr_arc arcToNext){
    assert(findArcMemberNoCompression(dec,arcToNext) == newRow->reducedMembers[reducedId].member);
    assert(SPQRarcIsValid(newRow->reducedMembers[reducedId].splitArc) && SPQRarcIsValid(arcToNext));
    SplitOrientation orientation;

    orientation.otherIsSource = newRow->reducedMembers[reducedId].otherIsSource;
    if(arcIsReversedNonRigid(dec,arcToNext) == arcIsReversedNonRigid(dec,newRow->reducedMembers[reducedId].splitArc)){
        orientation.headSplit = !newRow->reducedMembers[reducedId].splitHead;
    }else{
        orientation.headSplit = newRow->reducedMembers[reducedId].splitHead;
    }
    return orientation;
}

static SplitOrientation getRelativeOrientation(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                               reduced_member_id reducedId, spqr_arc arcToNext){
    switch(getMemberType(dec,newRow->reducedMembers[reducedId].member)){
        case SPQR_MEMBERTYPE_RIGID:
            return getRelativeOrientationRigid(dec,newRow,reducedId,arcToNext);
        case SPQR_MEMBERTYPE_PARALLEL:
            return getRelativeOrientationParallel(dec,newRow,reducedId,arcToNext);
        case SPQR_MEMBERTYPE_SERIES:
            return getRelativeOrientationSeries(dec,newRow,reducedId,arcToNext);
        default:
            assert(false);
    }
    SplitOrientation orientation;
    orientation.headSplit = false;
    orientation.otherIsSource = false;
    return orientation;
}
static void determineSplitTypeSeries(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                     reduced_member_id reducedId, spqr_arc marker,
                                     SplitOrientation previousOrientation){
    int numAdjacentMembers = newRow->reducedMembers[reducedId].numChildren-newRow->reducedMembers[reducedId].numPropagatedChildren;
    if(reducedMemberIsValid(newRow->reducedMembers[reducedId].parent) &&
       newRow->reducedMembers[newRow->reducedMembers[reducedId].parent].type != TYPE_PROPAGATED){
        ++numAdjacentMembers;
    }
    assert(numAdjacentMembers > 1);
    if(numAdjacentMembers > 2) {
        newRow->remainsNetwork = false;
        newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
        return;
    }
    cut_arc_id cutArc = newRow->reducedMembers[reducedId].firstCutArc;
    if(cutArcIsValid(cutArc)){
        spqr_arc arc = newRow->cutArcs[cutArc].arc;
        bool good = (((arcIsReversedNonRigid(dec,arc) == arcIsReversedNonRigid(dec,marker))
                     == newRow->cutArcs[cutArc].arcReversed ) == previousOrientation.headSplit)
                             == previousOrientation.otherIsSource;
        if(!good){
            newRow->remainsNetwork = false;
            newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
            return;
        }
        newRow->reducedMembers[reducedId].splitArc = marker;
        newRow->reducedMembers[reducedId].splitHead = previousOrientation.headSplit;
        newRow->reducedMembers[reducedId].otherIsSource = !previousOrientation.otherIsSource;
        newRow->reducedMembers[reducedId].type = TYPE_MERGED;
        return;
    }

    newRow->reducedMembers[reducedId].splitArc = marker;
    newRow->reducedMembers[reducedId].splitHead = previousOrientation.headSplit;
    newRow->reducedMembers[reducedId].otherIsSource = previousOrientation.otherIsSource;
    newRow->reducedMembers[reducedId].type = TYPE_MERGED;
}
static void determineSplitTypeParallel(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                     reduced_member_id reducedId, spqr_arc marker,
                                     SplitOrientation previousOrientation){
    MATRECRowReducedMember * redMember = &newRow->reducedMembers[reducedId];

    bool good = true;
    bool isReversed = true;
    int countedCutArcs = 0;
    for (cut_arc_id cutArc = redMember->firstCutArc; cutArcIsValid(cutArc);
         cutArc = newRow->cutArcs[cutArc].nextMember){
        spqr_arc arc = newRow->cutArcs[cutArc].arc;
        bool arcIsReversed = arcIsReversedNonRigid(dec,arc) != newRow->cutArcs[cutArc].arcReversed;
        if(countedCutArcs == 0){
            isReversed = arcIsReversed;
        }else if(arcIsReversed != isReversed){
            good = false;
            break;
        }
        ++countedCutArcs;
    }
    if(!good){
        redMember->type = TYPE_NOT_NETWORK;
        newRow->remainsNetwork = false;
        return;
    }
    if(countedCutArcs == 0){
        redMember->splitArc = marker;
        redMember->splitHead = previousOrientation.headSplit;
        redMember->otherIsSource = previousOrientation.otherIsSource;
        redMember->type = TYPE_MERGED;
        return;
    }
    bool isHeadSourceOrTailTarget = previousOrientation.headSplit == previousOrientation.otherIsSource;
    bool isOkay = isHeadSourceOrTailTarget == (isReversed == arcIsReversedNonRigid(dec,marker));
    if(!isOkay){
        redMember->type = TYPE_NOT_NETWORK;
        newRow->remainsNetwork = false;
        return;
    }
    redMember->splitArc = marker;
    redMember->splitHead = previousOrientation.headSplit;
    redMember->otherIsSource = previousOrientation.otherIsSource;
    redMember->type = TYPE_MERGED;
}
static void determineSplitTypeRigid(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                   reduced_member_id reducedId, spqr_member member, spqr_arc marker,
                                   SplitOrientation previousOrientation){
    assert(dec);
    assert(newRow);
    assert(getMemberType(dec,member) == SPQR_MEMBERTYPE_RIGID);

    Nodes nodes = rigidDetermineCandidateNodesFromAdjacentComponents(dec,newRow,reducedId);
    if(SPQRnodeIsInvalid(nodes.first) && SPQRnodeIsInvalid(nodes.second)){
        newRow->remainsNetwork = false;
        newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
        return;
    }
    if(SPQRnodeIsInvalid(nodes.first) && SPQRnodeIsValid(nodes.second)){
        nodes.first = nodes.second;
        nodes.second = SPQR_INVALID_NODE;
    }

    spqr_node markerHead = findArcHead(dec,marker);
    spqr_node markerTail = findArcTail(dec,marker);
    if(findArcSign(dec,marker).reversed){
        spqr_node temp = markerHead;
        markerHead = markerTail;
        markerTail = temp;
    }

    if(newRow->reducedMembers[reducedId].numCutArcs == 0){
        assert(SPQRnodeIsInvalid(nodes.second)); //There must be at least two adjacent components
        if(nodes.first != markerHead && nodes.first != markerTail){
            newRow->remainsNetwork = false;
            newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
            return;
        }
        bool reverse = previousOrientation.headSplit == (nodes.first == markerTail);
        newRow->reducedMembers[reducedId].splitNode = nodes.first;
        newRow->reducedMembers[reducedId].otherIsSource = previousOrientation.otherIsSource;
        newRow->reducedMembers[reducedId].willBeReversed = reverse;
        newRow->reducedMembers[reducedId].type = TYPE_MERGED;
        return;
    }
    if(!SPQRnodeIsValid(newRow->reducedMembers[reducedId].splitNode)){
        assert(newRow->reducedMembers[reducedId].numCutArcs > 0);//Checking for propagation only makes sense if there is at least one cut arc

        rigidFindStarNodes(dec,newRow,reducedId);
        if(!newRow->remainsNetwork){
            return;
        }

        if(SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode)){
            //not a star => attempt to find splittable nodes using articulation node algorithms
            //We save some work by telling the methods that only the marker nodes should be checked
            rigidGetSplittableArticulationPointsOnPath(dec,newRow,reducedId,nodes.first,nodes.second);
        }
        if(!newRow->remainsNetwork){
            return;
        }
        if(SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode)){
            newRow->remainsNetwork = false;
            newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
            return;
        }
    }
    spqr_node splitNode = determineAndColorSplitNode(dec,newRow,reducedId,nodes.first,nodes.second);

    if(SPQRnodeIsInvalid(splitNode)){
        newRow->remainsNetwork = false;
        newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
        return;
    }
    assert(splitNode == nodes.first || splitNode == nodes.second);
    assert(splitNode == markerHead || splitNode == markerTail);

    spqr_node otherNode = splitNode == markerHead ? markerTail : markerHead;
    bool headsMatch = previousOrientation.headSplit == (splitNode == markerHead);

    COLOR_STATUS otherColor = newRow->nodeColors[otherNode];
    assert(otherColor == COLOR_SOURCE || otherColor == COLOR_SINK);
    bool otherIsSource = otherColor == COLOR_SOURCE;

    bool good =  headsMatch == (previousOrientation.otherIsSource == otherIsSource);

    if(!good){
        newRow->remainsNetwork = false;
        newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
        return;
    }

    newRow->reducedMembers[reducedId].splitNode = splitNode;
    newRow->reducedMembers[reducedId].willBeReversed = !headsMatch;
    newRow->reducedMembers[reducedId].type = TYPE_MERGED;
}
static void determineSplitTypeNext(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                  reduced_member_id current, reduced_member_id next, bool currentIsParent){
    spqr_member member = newRow->reducedMembers[next].member;
    SPQRMemberType type = getMemberType(dec,member);
    spqr_arc nextArc = currentIsParent ? markerToParent(dec,member) : markerOfParent(dec,newRow->reducedMembers[current].member);
    spqr_arc currentArc = currentIsParent ? markerOfParent(dec,member) : markerToParent(dec,newRow->reducedMembers[current].member);
    SplitOrientation orientation = getRelativeOrientation(dec,newRow,current,currentArc);
    assert(type == SPQR_MEMBERTYPE_PARALLEL || type == SPQR_MEMBERTYPE_RIGID || type == SPQR_MEMBERTYPE_SERIES);
    switch(type){
        case SPQR_MEMBERTYPE_RIGID:{
            determineSplitTypeRigid(dec,newRow,next,member,nextArc,orientation);
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:{
            determineSplitTypeParallel(dec,newRow,next,nextArc,orientation);
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:{
            determineSplitTypeSeries(dec,newRow,next,nextArc,orientation);
            break;
        }
        default:
            assert(false);
            newRow->remainsNetwork = false;
    }
}
static void determineTypesChildrenNodes(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                        reduced_member_id parent, reduced_member_id node,
                                        reduced_member_id skipNode){
    if(node == skipNode || newRow->reducedMembers[node].type == TYPE_PROPAGATED){
        return;
    }
    //check merging
    determineSplitTypeNext(dec,newRow,parent,node,true);
    if(!newRow->remainsNetwork){
        return;
    }
    //merge all children
    for (int i = 0; i < newRow->reducedMembers[node].numChildren; ++i) {
        children_idx idx = newRow->reducedMembers[node].firstChild + i;
        reduced_member_id child = newRow->childrenStorage[idx];
        determineTypesChildrenNodes(dec,newRow,node,child,skipNode);
    }
}

static void determineMergeableTypes(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow, reduced_member_id root){
    assert(newRow->numReducedMembers <= newRow->memMergeTreeCallData);
    if(newRow->reducedMembers[root].numPropagatedChildren == newRow->reducedMembers[root].numChildren){
        //Determine single component;
        if(newRow->reducedMembers[root].type == TYPE_UNDETERMINED){
            spqr_member member = newRow->reducedMembers[root].member;
            switch(getMemberType(dec,member)){
                case SPQR_MEMBERTYPE_RIGID:
                    determineSingleRowRigidType(dec,newRow,root);
                    break;
                case SPQR_MEMBERTYPE_PARALLEL:
                    determineSingleParallelType(dec,newRow,root);
                    break;
                case SPQR_MEMBERTYPE_LOOP:
                case SPQR_MEMBERTYPE_SERIES:
                    determineSingleSeriesType(newRow,root);
                    break;
                default:
                    assert(false);
                    newRow->remainsNetwork = false;
            }
        }
        return;
    }

    //go to a leaf. We need to start in a leaf to avoid the ambiguity of choosing an orientation
    //in members which have no cut arcs; otherwise, we might choose the wrong one
    reduced_member_id leaf = root;
    while(newRow->reducedMembers[leaf].numChildren != newRow->reducedMembers[leaf].numPropagatedChildren){
        for(int i = 0; i < newRow->reducedMembers[leaf].numChildren;++i){
            children_idx idx = newRow->reducedMembers[leaf].firstChild + i;
            reduced_member_id child = newRow->childrenStorage[idx];
            if(newRow->reducedMembers[child].type != TYPE_PROPAGATED){
                leaf = child;
                break;
            }
        }
    }
    determineSplitTypeFirstLeaf(dec,newRow,leaf);

    if(!newRow->remainsNetwork){
        return;
    }
    reduced_member_id baseNode = leaf;
    reduced_member_id nextNode = newRow->reducedMembers[baseNode].parent;

    while(reducedMemberIsValid(nextNode)){
        //check this node
        determineSplitTypeNext(dec,newRow,baseNode,nextNode,false);
        if(!newRow->remainsNetwork){
            return;
        }
        //Add other nodes in the subtree
        for (int i = 0; i < newRow->reducedMembers[nextNode].numChildren; ++i) {
            children_idx idx = newRow->reducedMembers[nextNode].firstChild + i;
            reduced_member_id child = newRow->childrenStorage[idx];
            determineTypesChildrenNodes(dec,newRow,nextNode,child,baseNode);
            if(!newRow->remainsNetwork){
                return;
            }
        }

        //Move up one layer
        baseNode = nextNode;
        nextNode = newRow->reducedMembers[nextNode].parent;
    }
}

static void cleanUpRowMemberInformation(MATRECNetworkRowAddition * newRow){
    //This loop is at the end as memberInformation is also used to assign the cut arcs during propagation
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

static MATREC_ERROR rigidTransformArcIntoCycle(MATRECNetworkDecomposition *dec,
                                              const spqr_member member,
                                              const spqr_arc arc,
                                              const bool reverseArcDirection,
                                              NewRowInformation * const newRowInformation){
    //If a cycle already exists, just expand it with the new arc.
    spqr_member markerCycleMember = SPQR_INVALID_MEMBER;
    spqr_arc markerCycleArc = SPQR_INVALID_ARC;
    bool isParent = arc == markerToParent(dec,member);
    spqr_member adjacentMember = SPQR_INVALID_MEMBER;
    if (isParent) {
        adjacentMember = findMemberParent(dec, member);
        if (getMemberType(dec, adjacentMember) == SPQR_MEMBERTYPE_SERIES) {
            markerCycleMember = adjacentMember;
            markerCycleArc = markerOfParent(dec,member);
        }
    } else if (arcIsMarker(dec, arc)) {
        adjacentMember = findArcChildMember(dec, arc);
        if (getMemberType(dec, adjacentMember) == SPQR_MEMBERTYPE_SERIES) {
            markerCycleMember = adjacentMember;
            markerCycleArc = markerToParent(dec,adjacentMember);
        }
    }
    if (markerCycleMember != SPQR_INVALID_MEMBER) {
        newRowInformation->member = markerCycleMember;
        if(arcIsReversedNonRigid(dec,markerCycleArc)){
            newRowInformation->reversed = reverseArcDirection;
        }else{
            newRowInformation->reversed = !reverseArcDirection;
        }

        return MATREC_OKAY;
    }

    //Otherwise, we create a new cycle
    spqr_member newCycle;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &newCycle));
    //We would like to move the edge but unfortunately cannot do so without destroying the arc union-find datastructure.
    //Thus, we 'convert' the arc into a marker and add the new series

    spqr_arc duplicate = SPQR_INVALID_ARC;
    spqr_element element = arcGetElement(dec,arc);
    if(element != MARKER_COLUMN_ELEMENT && element != MARKER_ROW_ELEMENT){
        if(SPQRelementIsColumn(element)){
            MATREC_CALL(createColumnArc(dec,newCycle,&duplicate, SPQRelementToColumn(element),true));
        }else{
            MATREC_CALL(createRowArc(dec,newCycle,&duplicate, SPQRelementToRow(element),true));
        }
    }else if(isParent){
        //create parent marker
        MATREC_CALL(createParentMarker(dec, newCycle, arcIsTree(dec, arc), adjacentMember,
                                     markerOfParent(dec, member),&duplicate,true));
    }else{
        //create child marker
        MATREC_CALL(createChildMarker(dec,newCycle,adjacentMember,arcIsTree(dec,arc),&duplicate,true));
        dec->members[adjacentMember].parentMember = newCycle;
        dec->members[adjacentMember].markerOfParent = duplicate;
    }
        //Create the other marker edge
    spqr_arc cycleMarker = SPQR_INVALID_ARC;
    if(isParent){
        MATREC_CALL(createChildMarker(dec,newCycle,member,!arcIsTree(dec,arc),
                                    &cycleMarker,false));
    }else{
        MATREC_CALL(createParentMarker(dec,newCycle,!arcIsTree(dec,arc),
                                     member,arc,&cycleMarker,false));
    }
    //Change the existing edge to a marker
    if(isParent){
        assert(markerToParent(dec,member) == arc);
        dec->arcs[markerOfParent(dec,member)].childMember = newCycle;
        dec->members[member].parentMember = newCycle;
        dec->members[member].markerToParent = arc;
        dec->members[member].markerOfParent = cycleMarker;
        dec->arcs[arc].element =  arcIsTree(dec,arc) ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;;
        dec->arcs[arc].childMember = SPQR_INVALID_MEMBER;

    }else{
        dec->arcs[arc].element = arcIsTree(dec,arc) ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;
        dec->arcs[arc].childMember = newCycle;
    }
    newRowInformation->member = newCycle;
    newRowInformation->reversed = !reverseArcDirection;

    return MATREC_OKAY;
}

static MATREC_ERROR transformSingleRigid(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                       const reduced_member_id reducedMember,
                                       const spqr_member member,
                                       NewRowInformation * const newRowInformation){
    if(SPQRarcIsValid(newRow->reducedMembers[reducedMember].articulationArc)){
        spqr_arc arc = newRow->reducedMembers[reducedMember].articulationArc;
        //Cut arc is propagated into a cycle with new arc
        assert(newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcHeadNoCompression(dec,arc) ||
               newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcTailNoCompression(dec,arc));


        bool reversed;

        //TODO: check below...
        if(newRow->isArcCut[arc]){
            reversed = (newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcHead(dec,arc)) ==
                       newRow->reducedMembers[reducedMember].otherIsSource;
        }else{
            reversed = (newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcHead(dec,arc)) !=
                       newRow->reducedMembers[reducedMember].otherIsSource;
        }

        MATREC_CALL(rigidTransformArcIntoCycle(dec,member,newRow->reducedMembers[reducedMember].articulationArc,
                                             reversed, newRowInformation));

        return MATREC_OKAY;
    }
    //Single splittable node
    assert(SPQRnodeIsValid(newRow->reducedMembers[reducedMember].splitNode));

    spqr_node splitNode = newRow->reducedMembers[reducedMember].splitNode;
    if(newRow->reducedMembers[reducedMember].allHaveCommonNode){
        //Create a new node; move all cut arcs end of split node to it and add new arc between new node and split node
        spqr_node newNode = SPQR_INVALID_NODE;
        MATREC_CALL(createNode(dec, &newNode));

        cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
        do {
            spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
            spqr_node arcHead = findArcHead(dec, cutArc);
            if (arcHead == splitNode) {
                changeArcHead(dec, cutArc, arcHead, newNode);
            } else {
                changeArcTail(dec, cutArc, findArcTail(dec, cutArc), newNode);
            }

            cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
        } while (cutArcIsValid(cutArcIdx));

        newRowInformation->member = member;
        if(newRow->reducedMembers[reducedMember].otherIsSource){
            newRowInformation->head = newNode;
            newRowInformation->tail = splitNode;
        }else{
            newRowInformation->head = splitNode;
            newRowInformation->tail = newNode;
        }
        newRowInformation->representative = findArcSign(dec,newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc).representative;
        newRowInformation->reversed = false;

        return MATREC_OKAY;
    }
    //Articulation point was split (based on coloring)

    spqr_node newNode = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec, &newNode));

    spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
    spqr_arc iterArc = firstNodeArc;

    do{
        bool isCut = newRow->isArcCut[iterArc];
        spqr_node otherHead = findArcHead(dec, iterArc);
        spqr_node otherTail = findArcTail(dec, iterArc);
        spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
        bool isMoveColor = newRow->nodeColors[otherEnd] == COLOR_SOURCE;
        spqr_arc nextArc = getNextNodeArc(dec, iterArc, splitNode); //Need to do this before we modify the arc :)

        bool changeArcEnd = isCut == isMoveColor;
        if(changeArcEnd){
            if(otherHead == splitNode){
                changeArcHead(dec,iterArc,otherHead,newNode);
            }else{
                changeArcTail(dec,iterArc,otherTail,newNode);
            }
        }
        newRow->nodeColors[otherEnd] = UNCOLORED; //Clean up

        //Ugly hack to make sure we can iterate neighbourhood whilst changing arc ends.
        spqr_arc previousArc = iterArc;
        iterArc = nextArc;
        if(iterArc == firstNodeArc){
            break;
        }
        if(changeArcEnd && previousArc == firstNodeArc){
            firstNodeArc = iterArc;
        }
    }while(true);
    newRow->reducedMembers[reducedMember].coloredNode = SPQR_INVALID_NODE;

    newRowInformation->member = member;
    newRowInformation->head = newNode;
    newRowInformation->tail = splitNode;
    newRowInformation->representative = findArcSign(dec,newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc).representative;
    newRowInformation->reversed = false;

    return MATREC_OKAY;
}


static MATREC_ERROR splitParallelRowAddition(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                           const reduced_member_id reducedMember,
                                           const spqr_member member,
                                           NewRowInformation * newRowInfo){
    assert(newRow->reducedMembers[reducedMember].numCutArcs > 0);

    int numCutArcs = newRow->reducedMembers[reducedMember].numCutArcs;
    int numParallelArcs = getNumMemberArcs(dec,member);

    //TODO: also account for different tree arc orientation; if treeArc is reversed, everything is wrong.

    bool createCutParallel = numCutArcs > 1;
    bool convertOriginalParallel = (numCutArcs + 1) == numParallelArcs;

    //Do linear search to find non-marked arc
    spqr_arc treeArc = getFirstMemberArc(dec,member);
    do{
        if(arcIsTree(dec,treeArc)){
            break;
        }
        treeArc = getNextMemberArc(dec,treeArc);
    }while(treeArc != getFirstMemberArc(dec,member));
    assert(arcIsTree(dec,treeArc));

    bool treeReversed = arcIsReversedNonRigid(dec,treeArc);

    assert(!(!createCutParallel && convertOriginalParallel));//This can only happen if the parallel member is actually a loop, which means it is mislabeled
    if(createCutParallel) {
        if(convertOriginalParallel){
            spqr_member adjacentMember = SPQR_INVALID_MEMBER;
            spqr_arc adjacentArc = SPQR_INVALID_ARC;
            if(treeArc == markerToParent(dec,member)){
                adjacentMember = findMemberParent(dec,member);
                adjacentArc = markerOfParent(dec,member);
            }else if(arcIsMarker(dec,treeArc)){
                adjacentMember = findArcChildMember(dec,treeArc);
                adjacentArc = markerToParent(dec,adjacentMember);
                assert(markerOfParent(dec,adjacentMember) == treeArc);
            }
            cut_arc_id firstCut = newRow->reducedMembers[reducedMember].firstCutArc;
            bool firstReversed =  newRow->cutArcs[firstCut].arcReversed != arcIsReversedNonRigid(dec,newRow->cutArcs[firstCut].arc);

            if(SPQRmemberIsValid(adjacentMember) && getMemberType(dec, adjacentMember) == SPQR_MEMBERTYPE_SERIES){
                newRowInfo->member = adjacentMember;
                if(arcIsReversedNonRigid(dec,treeArc) == arcIsReversedNonRigid(dec,adjacentArc)){
                    newRowInfo->reversed = !firstReversed;
                }else{
                    newRowInfo->reversed = firstReversed;
                }
                return MATREC_OKAY;
            }
            spqr_member cutMember = SPQR_INVALID_MEMBER;
            MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

            cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
            assert(cutArcIsValid(cutArcIdx));
            bool parentCut = false;

            while(cutArcIsValid(cutArcIdx)){
                spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
                cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
                moveArcToNewMember(dec,cutArc,member,cutMember);
                if (cutArc == markerToParent(dec,member)){
                    parentCut = true;
                }
            }
            if(parentCut){
                MATREC_CALL(createMarkerPair(dec,cutMember,member,true,false,true));
            }else{
                MATREC_CALL(createMarkerPair(dec,member,cutMember,false,true,false));
            }
            changeLoopToSeries(dec,member);
            newRowInfo->member = member;
            if(treeReversed){
                newRowInfo->reversed = firstReversed == arcIsReversedNonRigid(dec,treeArc);
            }else{
                newRowInfo->reversed = firstReversed != arcIsReversedNonRigid(dec,treeArc);
            }
        }else{
            spqr_member cutMember = SPQR_INVALID_MEMBER;
            MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

            cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
            assert(cutArcIsValid(cutArcIdx));
            bool parentCut = false;

            while(cutArcIsValid(cutArcIdx)){
                spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
                cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
                moveArcToNewMember(dec,cutArc,member,cutMember);
                if (cutArc == markerToParent(dec,member)){
                    parentCut = true;
                }
            }
            spqr_member newSeries;
            MATREC_CALL(createMember(dec,SPQR_MEMBERTYPE_SERIES,&newSeries));
            if(parentCut){
                MATREC_CALL(createMarkerPair(dec,newSeries,member,true,false,false));
                MATREC_CALL(createMarkerPair(dec,cutMember,newSeries,true,false,true));
            }else{
                MATREC_CALL(createMarkerPair(dec,member,newSeries,false,false,false));
                MATREC_CALL(createMarkerPair(dec,newSeries,cutMember,false,true,false));
            }
            newRowInfo->member = newSeries;
            cut_arc_id firstCut = newRow->reducedMembers[reducedMember].firstCutArc;
            bool firstReversed =  newRow->cutArcs[firstCut].arcReversed != arcIsReversedNonRigid(dec,newRow->cutArcs[firstCut].arc);
            newRowInfo->reversed = firstReversed;
        }

        return MATREC_OKAY;
    }

    assert(!createCutParallel && !convertOriginalParallel);
    assert(numCutArcs == 1);
#ifndef NDEBUG
    spqr_arc arc = newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc;
    spqr_member adjacentMember = SPQR_INVALID_MEMBER;
    if(arc == markerToParent(dec,member)){
        adjacentMember = findMemberParent(dec,member);
    }else if(arcIsMarker(dec,arc)){
        adjacentMember = findArcChildMember(dec,arc);
    }
    if(SPQRmemberIsValid(adjacentMember)){
        assert(getMemberType(dec,adjacentMember) != SPQR_MEMBERTYPE_SERIES);
    }
#endif
    spqr_member newSeries;
    MATREC_CALL(createMember(dec,SPQR_MEMBERTYPE_SERIES,&newSeries));
    cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
    assert(cutArcIsValid(cutArcIdx));
    bool parentCut = false;

    while(cutArcIsValid(cutArcIdx)){
        spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
        cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
        moveArcToNewMember(dec,cutArc,member,newSeries);
        if (cutArc == markerToParent(dec,member)){
            parentCut = true;
        }
    }
    //TODO: fix parent/child orientations here
    if(parentCut){
        MATREC_CALL(createMarkerPair(dec,newSeries,member,true,true,false));
    }else{
        MATREC_CALL(createMarkerPair(dec,member,newSeries,false,false,true));
    }
    newRowInfo->member = newSeries;
    cut_arc_id cutArcId = newRow->reducedMembers[reducedMember].firstCutArc;
    newRowInfo->reversed = newRow->cutArcs[cutArcId].arcReversed == arcIsReversedNonRigid(dec,newRow->cutArcs[cutArcId].arc);

    return MATREC_OKAY;
}

static MATREC_ERROR transformSingleParallel(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                          const reduced_member_id reducedMember,
                                          const spqr_member member,
                                          NewRowInformation * info){
    MATREC_CALL(splitParallelRowAddition(dec,newRow,reducedMember,member,info));
    return MATREC_OKAY;
}

static MATREC_ERROR splitSeriesMergingRowAddition(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                                const reduced_member_id reducedMember,
                                                const spqr_member member,
                                                spqr_member * const mergingMember, bool * const isCut,
                                                spqr_arc * representativeEdge){
    assert(getNumMemberArcs(dec,member) >= 3);
    * isCut = newRow->reducedMembers[reducedMember].numCutArcs > 0;

    if(getNumMemberArcs(dec,member) == 3){
        spqr_arc splitArc = newRow->reducedMembers[reducedMember].splitArc;
        spqr_arc otherArc = SPQR_INVALID_ARC;
        for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
             i < newRow->reducedMembers[reducedMember].firstChild + newRow->reducedMembers[reducedMember].numChildren; ++i) {
            reduced_member_id child = newRow->childrenStorage[i];
            if(newRow->reducedMembers[child].type == TYPE_MERGED){
                spqr_arc testArc = markerOfParent(dec, findMember(dec,newRow->reducedMembers[child].member));
                if(testArc != splitArc){
                    otherArc = testArc;
                    break;
                }
            }
        }
        if(SPQRarcIsInvalid(otherArc)){
            reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
            assert(newRow->reducedMembers[parent].type == TYPE_MERGED || newRow->reducedMembers[parent].type == TYPE_PROPAGATED);
            assert(reducedMemberIsValid(parent) && newRow->reducedMembers[parent].type == TYPE_MERGED);
            otherArc = markerToParent(dec, member);

        }
        spqr_arc firstArc = getFirstMemberArc(dec,member);
        spqr_arc arc = firstArc;
        do{
            if(arc != splitArc && arc != otherArc){
                *representativeEdge = arc;
                break;
            }
            arc = getNextMemberArc(dec,arc);
        }while(arc != firstArc);
        *mergingMember = member;
        return MATREC_OKAY;
    }
    //Split off the relevant part of the series member
    spqr_member mergingSeries = SPQR_INVALID_MEMBER;
    MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &mergingSeries));
    //Move all marker arcs which point to another component in the reduced decomposition to the new member
    //This should be exactly 2, as with 3 the result is not network anymore
    //move all mergeable children and parent arcs to the mergingMember

    bool coTreeToMergingMember = false;
    bool parentToMergingMember = false;
    for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
         i < newRow->reducedMembers[reducedMember].firstChild + newRow->reducedMembers[reducedMember].numChildren; ++i) {
        reduced_member_id child = newRow->childrenStorage[i];
        assert(newRow->reducedMembers[child].type == TYPE_MERGED || newRow->reducedMembers[child].type == TYPE_PROPAGATED);
        if(newRow->reducedMembers[child].type == TYPE_MERGED){
            spqr_arc moveArc = markerOfParent(dec, findMember(dec,newRow->reducedMembers[child].member));
            moveArcToNewMember(dec,moveArc,member,mergingSeries);
            if(!arcIsTree(dec,moveArc)){
                coTreeToMergingMember = true;
            }
        }
    }

    reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
    assert(reducedMemberIsInvalid(parent) || (newRow->reducedMembers[parent].type == TYPE_MERGED || newRow->reducedMembers[parent].type == TYPE_PROPAGATED));

    if(reducedMemberIsValid(parent) &&
       newRow->reducedMembers[parent].type == TYPE_MERGED ){
        spqr_arc moveArc = markerToParent(dec, member);
        moveArcToNewMember(dec,moveArc,member,mergingSeries);
        parentToMergingMember = true;
        if(!arcIsTree(dec,moveArc)){
            coTreeToMergingMember = true;
        }
    }
    spqr_arc ignoreArc = SPQR_INVALID_ARC;
    if(parentToMergingMember){
        MATREC_CALL(createMarkerPairWithReferences(dec,mergingSeries,member,coTreeToMergingMember,true,false,representativeEdge,&ignoreArc));
    }else{
        MATREC_CALL(createMarkerPairWithReferences(dec,member,mergingSeries,!coTreeToMergingMember,false,true,&ignoreArc,representativeEdge));
    }

    *mergingMember = mergingSeries;
    assert(getNumMemberArcs(dec,mergingSeries) == 3 );
    assert(getNumMemberArcs(dec,member) >= 3);
    return MATREC_OKAY;
}

static MATREC_ERROR splitParallelMerging(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                       reduced_member_id reducedMember, spqr_member member,spqr_member * const pMergeMember,
                                       spqr_arc * const cutRepresentative){
    //When merging, we cannot have propagated members;
    assert(newRow->reducedMembers[reducedMember].numCutArcs < (getNumMemberArcs(dec,member)-1));

    int numMergeableAdjacent = newRow->reducedMembers[reducedMember].numChildren - newRow->reducedMembers[reducedMember].numPropagatedChildren;
    if(reducedMemberIsValid(newRow->reducedMembers[reducedMember].parent) &&
       newRow->reducedMembers[newRow->reducedMembers[reducedMember].parent].type == TYPE_MERGED){
        numMergeableAdjacent++;
    }

    int numCutArcs = newRow->reducedMembers[reducedMember].numCutArcs;
    //All arcs which are not in the mergeable decomposition or cut
    int numBaseSplitAwayArcs = getNumMemberArcs(dec,member) - numMergeableAdjacent - numCutArcs ;

    bool createCutParallel = numCutArcs > 1;
    bool keepOriginalParallel = numBaseSplitAwayArcs  <= 1;

    spqr_member cutMember = SPQR_INVALID_MEMBER;
    //The below cases can probably be aggregated in some way, but for now we first focus on the correct logic
    if(createCutParallel && keepOriginalParallel){
        bool parentCut = false;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

        cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
        assert(cutArcIsValid(cutArcIdx));

        while(cutArcIsValid(cutArcIdx)){
            spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
            cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
            moveArcToNewMember(dec,cutArc,member,cutMember);
            if (cutArc == markerToParent(dec,member)){
                parentCut = true;
            }
        }
        spqr_arc ignoreArc = SPQR_INVALID_ARC;
        if(parentCut){
            MATREC_CALL(createMarkerPairWithReferences(dec,cutMember,member,true,false,false,&ignoreArc,cutRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,cutMember,false,false,false,cutRepresentative,&ignoreArc));
        }

        *pMergeMember = member;
    }else if(createCutParallel){
        assert(!keepOriginalParallel);

        bool parentCut = false;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

        cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
        assert(cutArcIsValid(cutArcIdx));

        while(cutArcIsValid(cutArcIdx)){
            spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
            cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
            moveArcToNewMember(dec,cutArc,member,cutMember);
            if (cutArc == markerToParent(dec,member)){
                parentCut = true;
            }
        }
        spqr_arc ignoreArc = SPQR_INVALID_ARC;
        if(parentCut){
            MATREC_CALL(createMarkerPairWithReferences(dec,cutMember,member,true,false,false,&ignoreArc,cutRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,cutMember,false,false,false,cutRepresentative,&ignoreArc));
        }


        spqr_arc noCutRepresentative = SPQR_INVALID_ARC;
        spqr_member mergingMember = member;
        bool parentToMergingMember = false;
        bool treeToMergingMember = false;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &mergingMember));
        //move all mergeable children and parent arcs to the mergingMember
        for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
             i < newRow->reducedMembers[reducedMember].firstChild +
                 newRow->reducedMembers[reducedMember].numChildren; ++i) {
            reduced_member_id child = newRow->childrenStorage[i];
            assert(newRow->reducedMembers[child].type == TYPE_MERGED || newRow->reducedMembers[child].type == TYPE_PROPAGATED);
            if (newRow->reducedMembers[child].type == TYPE_MERGED) {
                spqr_arc moveArc = markerOfParent(dec, findMember(dec,newRow->reducedMembers[child].member));
                moveArcToNewMember(dec, moveArc, member, mergingMember);
                if (arcIsTree(dec, moveArc)) {
                    treeToMergingMember = true;
                }
            }
        }
        reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
        if (reducedMemberIsValid(parent) &&
            newRow->reducedMembers[parent].type == TYPE_MERGED) {
            spqr_arc moveArc = markerToParent(dec, member);
            moveArcToNewMember(dec, moveArc, member, mergingMember);
            parentToMergingMember = true;
            if (arcIsTree(dec, moveArc)) {
                treeToMergingMember = true;
            }
        }
        //If there is only one cut arc, we also move it.
        if (SPQRarcIsValid(*cutRepresentative)) {
            if (*cutRepresentative == markerToParent(dec, member)) {
                parentToMergingMember = true;
            }
            moveArcToNewMember(dec, *cutRepresentative, member, mergingMember);
        }
        spqr_arc ignoreArgument = SPQR_INVALID_ARC;
        if(parentToMergingMember || parentCut){
            MATREC_CALL(createMarkerPairWithReferences(dec,mergingMember,member,!treeToMergingMember,false,false,&ignoreArgument,&noCutRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,mergingMember,treeToMergingMember,false,false,&noCutRepresentative,&ignoreArgument));
        }
        *pMergeMember = mergingMember;
    }else if(keepOriginalParallel){
        assert(!createCutParallel);
        if(cutArcIsValid( newRow->reducedMembers[reducedMember].firstCutArc)){
            *cutRepresentative = newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc;
        }
        *pMergeMember = member;

    }else{
        assert(!keepOriginalParallel && ! createCutParallel);

        if(cutArcIsValid( newRow->reducedMembers[reducedMember].firstCutArc)){
            *cutRepresentative = newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc;
        }

        spqr_arc noCutRepresentative = SPQR_INVALID_ARC;
        spqr_member mergingMember = member;
        bool parentToMergingMember = false;
        bool treeToMergingMember = false;
        MATREC_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &mergingMember));
        //move all mergeable children and parent arcs to the mergingMember
        for (children_idx i = newRow->reducedMembers[reducedMember].firstChild;
             i < newRow->reducedMembers[reducedMember].firstChild +
                 newRow->reducedMembers[reducedMember].numChildren; ++i) {
            reduced_member_id child = newRow->childrenStorage[i];
            assert(newRow->reducedMembers[child].type == TYPE_MERGED || newRow->reducedMembers[child].type == TYPE_PROPAGATED);
            if (newRow->reducedMembers[child].type == TYPE_MERGED) {
                spqr_arc moveArc = markerOfParent(dec, findMember(dec,newRow->reducedMembers[child].member));
                moveArcToNewMember(dec, moveArc, member, mergingMember);
                if (arcIsTree(dec, moveArc)) {
                    treeToMergingMember = true;
                }
            }
        }
        reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
        if (reducedMemberIsValid(parent) &&
            newRow->reducedMembers[parent].type == TYPE_MERGED) {
            spqr_arc moveArc = markerToParent(dec, member);
            moveArcToNewMember(dec, moveArc, member, mergingMember);
            parentToMergingMember = true;
            if (arcIsTree(dec, moveArc)) {
                treeToMergingMember = true;
            }
        }
        //If there is only one cut arc, we also move it.
        if (SPQRarcIsValid(*cutRepresentative)) {
            if (*cutRepresentative == markerToParent(dec, member)) {
                parentToMergingMember = true;
            }
            moveArcToNewMember(dec, *cutRepresentative, member, mergingMember);
        }
        spqr_arc ignoreArgument = SPQR_INVALID_ARC;
        if(parentToMergingMember){
            MATREC_CALL(createMarkerPairWithReferences(dec,mergingMember,member,!treeToMergingMember,false,false,&ignoreArgument,&noCutRepresentative));
        }else{
            MATREC_CALL(createMarkerPairWithReferences(dec,member,mergingMember,treeToMergingMember,false,false,&noCutRepresentative,&ignoreArgument));
        }
        *pMergeMember = mergingMember;
    }
    return MATREC_OKAY;
}
static MATREC_ERROR splitFirstLeaf(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                 reduced_member_id leaf, NewRowInformation * const newRowInformation){
    assert(cutArcIsValid(newRow->reducedMembers[leaf].firstCutArc));
    assert(newRow->reducedMembers[leaf].numCutArcs > 0);
    spqr_member member = newRow->reducedMembers[leaf].member;
    if(getMemberType(dec,member) == SPQR_MEMBERTYPE_PARALLEL){
        spqr_member mergeMember = SPQR_INVALID_MEMBER;
        spqr_arc cutRepresentative = SPQR_INVALID_ARC;
        MATREC_CALL(splitParallelMerging(dec,newRow,leaf,member,&mergeMember,&cutRepresentative));
        newRow->reducedMembers[leaf].member = mergeMember;

        spqr_node firstNode = SPQR_INVALID_NODE;
        spqr_node secondNode = SPQR_INVALID_NODE;
        spqr_node thirdNode = SPQR_INVALID_NODE;
        MATREC_CALL(createNode(dec,&firstNode));
        MATREC_CALL(createNode(dec,&secondNode));
        MATREC_CALL(createNode(dec,&thirdNode));

        spqr_arc splitArc = newRow->reducedMembers[leaf].splitArc;
        assert(findArcMemberNoCompression(dec,splitArc) == mergeMember);
        bool splitArcReversed = arcIsReversedNonRigid(dec,splitArc);
        bool splitHead = newRow->reducedMembers[leaf].splitHead;
        spqr_node splitArcHead = splitArcReversed ? secondNode : firstNode;
        spqr_node splitArcTail = splitArcReversed ? firstNode : secondNode;
        spqr_node splitNode = splitHead ? splitArcHead : splitArcTail;
        spqr_node otherNode = splitHead ? splitArcTail : splitArcHead;

        spqr_arc first_arc = getFirstMemberArc(dec, mergeMember);
        spqr_arc arc = first_arc;

        do {
            if(arc != cutRepresentative){
                if(arcIsReversedNonRigid(dec,arc)){
                    setArcHeadAndTail(dec,arc,secondNode,firstNode);
                }else{
                    setArcHeadAndTail(dec,arc,firstNode,secondNode);
                }
            }else{
                if ((arcIsReversedNonRigid(dec, arc) == splitArcReversed) == splitHead) {
                    setArcHeadAndTail(dec, arc, thirdNode, otherNode);
                } else {
                    setArcHeadAndTail(dec, arc, otherNode, thirdNode);
                }
            }
            arcSetReversed(dec,arc,false);
            if(arc == splitArc){
                arcSetRepresentative(dec,arc,SPQR_INVALID_ARC);
            }else{
                arcSetRepresentative(dec,arc,splitArc);
            }
            arc = getNextMemberArc(dec, arc);
        } while (arc != first_arc);


        updateMemberType(dec, mergeMember, SPQR_MEMBERTYPE_RIGID);
        newRowInformation->member = mergeMember;
        if(newRow->reducedMembers[leaf].otherIsSource){
            newRowInformation->head = thirdNode;
            newRowInformation->tail = splitNode;
        }else{
            newRowInformation->head = splitNode;
            newRowInformation->tail = thirdNode;
        }

        newRowInformation->reversed = false;
        newRowInformation->representative = splitArc;

        return MATREC_OKAY;
    }
    assert(getMemberType(dec,member) == SPQR_MEMBERTYPE_RIGID);

    spqr_node newNode = SPQR_INVALID_NODE; //Sink node
    MATREC_CALL(createNode(dec,&newNode));

    spqr_node splitNode = newRow->reducedMembers[leaf].splitNode;

    spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
    spqr_arc iterArc = firstNodeArc;
    do {
        spqr_node otherHead = findArcHead(dec, iterArc);
        spqr_node otherTail = findArcTail(dec, iterArc);
        spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
        spqr_arc nextArc = getNextNodeArc(dec, iterArc, splitNode); //Need to do this before we modify the arc

        bool isCut = newRow->isArcCut[iterArc];
        bool isMoveColor = newRow->nodeColors[otherEnd] == COLOR_SOURCE;
        bool changeArcEnd = isCut == isMoveColor;
        if (changeArcEnd) {
            if (otherHead == splitNode) {
                changeArcHead(dec, iterArc, otherHead, newNode);
            } else {
                changeArcTail(dec, iterArc, otherTail, newNode);
            }
        }
        //Ugly hack to make sure we can iterate neighbourhood whilst changing arc ends.
        newRow->nodeColors[otherEnd] = UNCOLORED;
        spqr_arc previousArc = iterArc;
        iterArc = nextArc;
        if (iterArc == firstNodeArc) {
            break;
        }
        if (changeArcEnd && previousArc == firstNodeArc) {
            firstNodeArc = iterArc;
        }
    } while (true);
    newRowInformation->head = newNode;
    newRowInformation->tail = splitNode;
    newRowInformation->member = member;
    newRowInformation->reversed = false;
    newRowInformation->representative = findArcSign(dec,iterArc).representative;

    return MATREC_OKAY;
}

static MATREC_ERROR mergeSplitMemberIntoParent(MATRECNetworkDecomposition *dec,
                                             MATRECNetworkRowAddition * newRow,
                                             spqr_member member,
                                             spqr_member parent,
                                             spqr_arc parentToChild,
                                             spqr_arc childToParent,
                                             bool headToHead,
                                             spqr_node parentNode,
                                             spqr_node childNode,
                                             spqr_member * mergedMember,
                                             spqr_node * arcNodeOne,
                                             spqr_node * arcNodeTwo,
                                             spqr_node * thirdNode) {
    assert(dec);
    assert(SPQRmemberIsValid(member));
    assert(memberIsRepresentative(dec, member));
    assert(SPQRmemberIsValid(parent));
    assert(memberIsRepresentative(dec, parent));
    assert(findMemberParentNoCompression(dec, member) == parent);
    assert(markerOfParent(dec, member) == parentToChild);
    assert(markerToParent(dec, member) == childToParent);

    removeArcFromMemberArcList(dec, parentToChild, parent);
    removeArcFromMemberArcList(dec, childToParent, member);

    spqr_node parentArcNodes[2] = {findEffectiveArcTail(dec, parentToChild), findEffectiveArcHead(dec, parentToChild)};
    spqr_node childArcNodes[2] = {findEffectiveArcTail(dec, childToParent), findEffectiveArcHead(dec, childToParent)};

    clearArcHeadAndTail(dec, parentToChild);
    clearArcHeadAndTail(dec, childToParent);

    spqr_node first = childArcNodes[headToHead ? 0 : 1];
    spqr_node second = childArcNodes[headToHead ? 1 : 0];
    {
        spqr_node newNode = mergeNodes(dec, parentArcNodes[0], first);
        spqr_node toRemoveFrom = newNode == first ? parentArcNodes[0] : first;
        mergeNodeArcList(dec, newNode, toRemoveFrom);
        *arcNodeOne = newNode;
        newRow->nodeColors[toRemoveFrom] = UNCOLORED;
    }
    {
        spqr_node newNode = mergeNodes(dec, parentArcNodes[1], second);
        spqr_node toRemoveFrom = newNode == second ? parentArcNodes[1] : second;
        mergeNodeArcList(dec, newNode, toRemoveFrom);
        *arcNodeTwo = newNode;
        newRow->nodeColors[toRemoveFrom] = UNCOLORED;
    }
    {
        spqr_node newNode = mergeNodes(dec, parentNode, childNode);
        spqr_node toRemoveFrom = newNode == childNode ? parentNode : childNode;
        mergeNodeArcList(dec,newNode,toRemoveFrom);
        *thirdNode = newNode;
        newRow->nodeColors[toRemoveFrom] = UNCOLORED;
    }

    spqr_member newMember = mergeMembers(dec, member, parent);
    spqr_member toRemoveFrom = newMember == member ? parent : member;
    mergeMemberArcList(dec, newMember, toRemoveFrom);
    if (toRemoveFrom == parent) {
        updateMemberParentInformation(dec, newMember, toRemoveFrom);
    }
    updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);
    *mergedMember = newMember;
    return MATREC_OKAY;
}
static MATREC_ERROR splitAndMergeSeries(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                       reduced_member_id smallMember,
                                      bool largeIsParent, NewRowInformation * const newRowInformation,
                                      spqr_member member){
    bool isCut = false;
    spqr_member mergingMember = SPQR_INVALID_MEMBER;
    spqr_arc nonVirtualArc = SPQR_INVALID_ARC;
    MATREC_CALL(splitSeriesMergingRowAddition(dec,newRow,smallMember,member,&mergingMember,&isCut,&nonVirtualArc));
    assert(getNumMemberArcs(dec,mergingMember) == 3);

    //create the split series. There's two possible configurations, based on whether it contains a cut edge or not
    spqr_node a = SPQR_INVALID_NODE;
    spqr_node b = SPQR_INVALID_NODE;
    spqr_node c = SPQR_INVALID_NODE;
    spqr_node d = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec,&a));
    MATREC_CALL(createNode(dec,&b));
    MATREC_CALL(createNode(dec,&c));
    MATREC_CALL(createNode(dec,&d));

    spqr_arc splitArc = newRow->reducedMembers[smallMember].splitArc;

    {
        bool splitHead = newRow->reducedMembers[smallMember].splitHead;
        spqr_arc firstArc = getFirstMemberArc(dec,mergingMember);
        spqr_arc arc = firstArc;
        bool splitReversed = arcIsReversedNonRigid(dec,splitArc);
        do{
            if(arc == splitArc){
                if(splitHead){
                    setArcHeadAndTail(dec,splitArc,b,a);
                }else{
                    setArcHeadAndTail(dec,splitArc,a,b);
                }
            }else if (arc == nonVirtualArc){
                if ((arcIsReversedNonRigid(dec, arc) == splitReversed) == splitHead) {
                    setArcHeadAndTail(dec, arc, a, d);
                } else {
                    setArcHeadAndTail(dec, arc, d, a);
                }
            }else{
                spqr_node otherNode = cutArcIsValid(newRow->reducedMembers[smallMember].firstCutArc) ? c : b;
                if((arcIsReversedNonRigid(dec,arc) == splitReversed) == splitHead){
                    setArcHeadAndTail(dec, arc, d, otherNode);
                }else{
                    setArcHeadAndTail(dec, arc, otherNode, d);
                }
            }
            arcSetReversed(dec,arc,false);
            arcSetRepresentative(dec,arc,splitArc);
            arc = getNextMemberArc(dec,arc);
        }while(arc != firstArc);
        arcSetRepresentative(dec,splitArc,SPQR_INVALID_ARC);
    }

    spqr_member otherMember = newRowInformation->member;
    spqr_arc otherMarker = largeIsParent ? markerOfParent(dec,mergingMember) : markerToParent(dec,otherMember);

    assert(nodeIsRepresentative(dec,newRowInformation->tail));
    assert(nodeIsRepresentative(dec,newRowInformation->head));
    spqr_node splitNode = newRow->reducedMembers[smallMember].splitHead ? findEffectiveArcHead(dec,otherMarker) :
            findEffectiveArcTail(dec,otherMarker);

    spqr_node otherNode = splitNode == newRowInformation->head ? newRowInformation->tail : newRowInformation->head;
    assert(splitNode == newRowInformation->head || splitNode == newRowInformation->tail);
    newRowInformation->representative =  mergeArcSigns(dec,newRowInformation->representative,splitArc,false);

    spqr_member mergedMember = SPQR_INVALID_MEMBER;
    spqr_node arcNodeOne;
    spqr_node arcNodeTwo;
    spqr_node thirdNode;
    if(largeIsParent){
        MATREC_CALL(mergeSplitMemberIntoParent(dec,newRow,mergingMember,otherMember,otherMarker,splitArc,true,
                                   otherNode,c,&mergedMember,
                                   &arcNodeOne,
                                   &arcNodeTwo,
                                   &thirdNode));
    }else{
        MATREC_CALL(mergeSplitMemberIntoParent(dec,newRow,otherMember,mergingMember, splitArc, otherMarker,true,
                                             c,otherNode,&mergedMember,
                                             &arcNodeOne,
                                             &arcNodeTwo,
                                             &thirdNode));
    }
    newRow->reducedMembers[smallMember].member = mergedMember;

    newRowInformation->member = mergedMember;
    bool splitIsReferenceHead = newRow->reducedMembers[smallMember].splitHead;
    bool splitIsNewRowHead = splitNode == newRowInformation->head;
    if(!splitIsReferenceHead && !splitIsNewRowHead){
        newRowInformation->head = thirdNode;
        newRowInformation->tail = arcNodeOne;
    }else if (!splitIsReferenceHead && splitIsNewRowHead){
        newRowInformation->head = arcNodeOne;
        newRowInformation->tail = thirdNode;
    }else if(splitIsReferenceHead && !splitIsNewRowHead){
        newRowInformation->head = thirdNode;
        newRowInformation->tail = arcNodeTwo;
    }else if (splitIsReferenceHead && splitIsNewRowHead){
        newRowInformation->head = arcNodeTwo;
        newRowInformation->tail = thirdNode;
    }

    return MATREC_OKAY;
}

static MATREC_ERROR splitAndMergeParallel(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                          reduced_member_id smallMember,
                                      bool largeIsParent, NewRowInformation * const newRowInformation,
                                      spqr_member member){
    spqr_member mergeMember = SPQR_INVALID_MEMBER;
    spqr_arc cutRepresentative = SPQR_INVALID_ARC;
    MATREC_CALL(splitParallelMerging(dec,newRow,smallMember,member,&mergeMember,&cutRepresentative));
    newRow->reducedMembers[smallMember].member = mergeMember;

    spqr_node firstNode = SPQR_INVALID_NODE;
    spqr_node secondNode = SPQR_INVALID_NODE;
    spqr_node thirdNode = SPQR_INVALID_NODE;
    MATREC_CALL(createNode(dec,&firstNode));
    MATREC_CALL(createNode(dec,&secondNode));
    MATREC_CALL(createNode(dec,&thirdNode));

    spqr_arc splitArc = newRow->reducedMembers[smallMember].splitArc;
    assert(findArcMemberNoCompression(dec,splitArc) == mergeMember);
    bool splitArcReversed = arcIsReversedNonRigid(dec,splitArc);
    bool splitHead = newRow->reducedMembers[smallMember].splitHead;

    spqr_node splitArcHead = splitArcReversed ? secondNode : firstNode;
    spqr_node splitArcTail = splitArcReversed ? firstNode : secondNode;
    spqr_node otherNode = splitHead ? splitArcTail : splitArcHead;

    spqr_arc first_arc = getFirstMemberArc(dec, mergeMember);
    spqr_arc arc = first_arc;

    do {
        if(arc != cutRepresentative){
            if(arcIsReversedNonRigid(dec,arc)){
                setArcHeadAndTail(dec,arc,secondNode,firstNode);
            }else{
                setArcHeadAndTail(dec,arc,firstNode,secondNode);
            }
        }else{
            if ((arcIsReversedNonRigid(dec, arc) == splitArcReversed) == splitHead) {
                setArcHeadAndTail(dec, arc, thirdNode, otherNode);
            } else {
                setArcHeadAndTail(dec, arc, otherNode, thirdNode);
            }
        }
        arcSetReversed(dec,arc,false);
        arcSetRepresentative(dec,arc,splitArc);
        arc = getNextMemberArc(dec, arc);
    } while (arc != first_arc);
    arcSetRepresentative(dec,splitArc,SPQR_INVALID_ARC);

    spqr_member otherMember = newRowInformation->member;
    spqr_arc otherMarker = largeIsParent ? markerOfParent(dec,mergeMember) : markerToParent(dec,otherMember);

    assert(nodeIsRepresentative(dec,newRowInformation->tail));
    assert(nodeIsRepresentative(dec,newRowInformation->head));
    spqr_node largeSplitNode = newRow->reducedMembers[smallMember].splitHead ? findEffectiveArcHead(dec,otherMarker) :
                          findEffectiveArcTail(dec,otherMarker);

    spqr_node largeOtherNode = largeSplitNode == newRowInformation->head ? newRowInformation->tail : newRowInformation->head;
    assert(largeSplitNode == newRowInformation->head || largeSplitNode == newRowInformation->tail);

    newRowInformation->representative = mergeArcSigns(dec,newRowInformation->representative,splitArc,false);

    spqr_member mergedMember = SPQR_INVALID_MEMBER;
    spqr_node arcNodeOne;
    spqr_node arcNodeTwo;
    spqr_node mergeNodeThree;
    if(largeIsParent){
        MATREC_CALL(mergeSplitMemberIntoParent(dec,newRow,mergeMember,otherMember,otherMarker,splitArc,true,
                                             largeOtherNode,thirdNode,&mergedMember,
                                             &arcNodeOne,
                                             &arcNodeTwo,
                                             &mergeNodeThree));
    }else{
        MATREC_CALL(mergeSplitMemberIntoParent(dec,newRow,otherMember,mergeMember, splitArc, otherMarker,true,
                                             thirdNode,largeOtherNode,&mergedMember,
                                             &arcNodeOne,
                                             &arcNodeTwo,
                                             &mergeNodeThree));
    }


    newRowInformation->member = mergedMember;

    bool splitIsReferenceHead = newRow->reducedMembers[smallMember].splitHead;
    bool splitIsNewRowHead = largeSplitNode == newRowInformation->head;
    if(!splitIsReferenceHead && !splitIsNewRowHead){
        newRowInformation->head = mergeNodeThree;
        newRowInformation->tail = arcNodeOne;
    }else if (!splitIsReferenceHead && splitIsNewRowHead){
        newRowInformation->head = arcNodeOne;
        newRowInformation->tail = mergeNodeThree;
    }else if(splitIsReferenceHead && !splitIsNewRowHead){
        newRowInformation->head = mergeNodeThree;
        newRowInformation->tail = arcNodeTwo;
    }else if (splitIsReferenceHead && splitIsNewRowHead){
        newRowInformation->head = arcNodeTwo;
        newRowInformation->tail = mergeNodeThree;
    }

    return MATREC_OKAY;
}
static MATREC_ERROR splitAndMergeRigid(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                       reduced_member_id smallMember,
                                     bool largeIsParent, NewRowInformation * const newRowInformation,
                                     spqr_member member){

    spqr_node newNode = SPQR_INVALID_NODE; //Sink node
    MATREC_CALL(createNode(dec,&newNode));

    spqr_member smallMemberMember = member;
    spqr_member largeMemberMember = newRowInformation->member;

    spqr_arc smallMarker = largeIsParent ? markerToParent(dec,smallMemberMember) : markerOfParent(dec,largeMemberMember);
    spqr_arc largeMarker = largeIsParent ? markerOfParent(dec,smallMemberMember) : markerToParent(dec,largeMemberMember);

    spqr_node splitNode = newRow->reducedMembers[smallMember].splitNode;
    spqr_node smallOtherNode = newNode;

    if(newRow->reducedMembers[smallMember].numCutArcs != 0) {
        spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
        spqr_arc iterArc = firstNodeArc;
        do {
            spqr_node otherHead = findArcHead(dec, iterArc);
            spqr_node otherTail = findArcTail(dec, iterArc);
            spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
            spqr_arc nextArc = getNextNodeArc(dec, iterArc, splitNode); //Need to do this before we modify the arc

            bool isCut = newRow->isArcCut[iterArc];
            bool isMoveColor = newRow->nodeColors[otherEnd] == COLOR_SOURCE;
            bool changeArcEnd = isCut == isMoveColor;
            if (changeArcEnd) {
                if (otherHead == splitNode) {
                    changeArcHead(dec, iterArc, otherHead, newNode);
                } else {
                    changeArcTail(dec, iterArc, otherTail, newNode);
                }
                if(iterArc == smallMarker){
                    smallOtherNode = splitNode;
                }
            }
            newRow->nodeColors[otherEnd] = UNCOLORED;
            //Ugly hack to make sure we can iterate neighbourhood whilst changing arc ends.
            spqr_arc previousArc = iterArc;
            iterArc = nextArc;
            if (iterArc == firstNodeArc) {
                break;
            }
            if (changeArcEnd && previousArc == firstNodeArc) {
                firstNodeArc = iterArc;
            }
        } while (true);
    }

    spqr_arc representative = findArcSign(dec,smallMarker).representative;

    newRowInformation->representative = mergeArcSigns(dec,newRowInformation->representative,representative,
                                                      newRow->reducedMembers[smallMember].willBeReversed);

    spqr_node largeMarkerHead = findArcHead(dec,largeMarker);
    spqr_node largeMarkerTail = findArcTail(dec,largeMarker);
    if(findArcSign(dec,largeMarker).reversed){
        spqr_node temp = largeMarkerHead;
        largeMarkerHead = largeMarkerTail;
        largeMarkerTail = temp;
    }
    assert(newRowInformation->head == largeMarkerHead || newRowInformation->head == largeMarkerTail ||
           newRowInformation->tail == largeMarkerHead || newRowInformation->tail == largeMarkerTail);
    spqr_node largeOtherNode = (newRowInformation->head == largeMarkerHead || newRowInformation->head == largeMarkerTail)
                               ? newRowInformation->tail : newRowInformation->head;

    spqr_member mergedMember = SPQR_INVALID_MEMBER;
    spqr_node arcNodeOne;
    spqr_node arcNodeTwo;
    spqr_node mergeNodeThree;
    if(largeIsParent){
        MATREC_CALL(mergeSplitMemberIntoParent(dec,newRow,smallMemberMember,largeMemberMember,largeMarker,smallMarker,true,
                                             largeOtherNode,smallOtherNode,&mergedMember,
                                             &arcNodeOne,
                                             &arcNodeTwo,
                                             &mergeNodeThree));
    }else{
        MATREC_CALL(mergeSplitMemberIntoParent(dec,newRow,largeMemberMember,smallMemberMember,smallMarker,largeMarker,true,
                                             smallOtherNode,largeOtherNode,&mergedMember,
                                             &arcNodeOne,
                                             &arcNodeTwo,
                                             &mergeNodeThree));
    }
    newRowInformation->member = mergedMember;

    bool otherIsHead = largeOtherNode == newRowInformation->head;
    bool adjacentToMarkerHead = (newRowInformation->tail == largeMarkerHead || newRowInformation->head == largeMarkerHead);
    if(adjacentToMarkerHead){
        if(otherIsHead){
            newRowInformation->head = mergeNodeThree;
            newRowInformation->tail = arcNodeTwo;
        }else{
            newRowInformation->head = arcNodeTwo;
            newRowInformation->tail = mergeNodeThree;
        }
    }else{
        if(otherIsHead){
            newRowInformation->head = mergeNodeThree;
            newRowInformation->tail = arcNodeOne;
        }else{
            newRowInformation->head = arcNodeOne;
            newRowInformation->tail = mergeNodeThree;
        }
    }

    return MATREC_OKAY;
}
static MATREC_ERROR splitAndMerge(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                 reduced_member_id smallMember,
                                bool largeIsParent, NewRowInformation * const newRowInformation){
    spqr_member member = newRow->reducedMembers[smallMember].member;
    switch(getMemberType(dec,member)){
        case SPQR_MEMBERTYPE_RIGID:{
            MATREC_CALL(splitAndMergeRigid(dec,newRow,smallMember,largeIsParent,newRowInformation,member));
            break;
        }
        case SPQR_MEMBERTYPE_PARALLEL:{
            MATREC_CALL(splitAndMergeParallel(dec,newRow,smallMember,largeIsParent,newRowInformation,member));
            break;
        }
        case SPQR_MEMBERTYPE_SERIES:{
            MATREC_CALL(splitAndMergeSeries(dec,newRow,smallMember,largeIsParent,newRowInformation,member));
            break;
        }
        default:
            assert(false);
            newRow->remainsNetwork = false;
    }
    return MATREC_OKAY;
}
static MATREC_ERROR mergeChildrenNodes(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
                                      reduced_member_id node,
                                     reduced_member_id skipNode, NewRowInformation * const newRowInformation
                                     ){
    if(node == skipNode || newRow->reducedMembers[node].type == TYPE_PROPAGATED){
        return MATREC_OKAY;
    }
    //check merging
    MATREC_CALL(splitAndMerge(dec,newRow,node,true,newRowInformation));

    //merge all children
    for (int i = 0; i < newRow->reducedMembers[node].numChildren; ++i) {
        children_idx idx = newRow->reducedMembers[node].firstChild + i;
        reduced_member_id child = newRow->childrenStorage[idx];
        MATREC_CALL(mergeChildrenNodes(dec,newRow,child,skipNode,newRowInformation));
    }
    return MATREC_OKAY;
}
static MATREC_ERROR mergeTree(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow, reduced_member_id root,
                            NewRowInformation * const newRowInformation){

    //We use the same ordering as when finding
    //go to a leaf. We need to start in a leaf to avoid the ambiguity of choosing an orientation
    //in members which have no cut arcs; otherwise, we might choose the wrong one
    reduced_member_id leaf = root;
    while(newRow->reducedMembers[leaf].numChildren != newRow->reducedMembers[leaf].numPropagatedChildren){
        for(int i = 0; i < newRow->reducedMembers[leaf].numChildren;++i){
            children_idx idx = newRow->reducedMembers[leaf].firstChild + i;
            reduced_member_id child = newRow->childrenStorage[idx];
            if(newRow->reducedMembers[child].type != TYPE_PROPAGATED){
                leaf = child;
                break;
            }
        }
    }
    MATREC_CALL(splitFirstLeaf(dec,newRow,leaf,newRowInformation));

    reduced_member_id baseNode = leaf;
    reduced_member_id nextNode = newRow->reducedMembers[baseNode].parent;

    while(reducedMemberIsValid(nextNode)){
        //check this node
        MATREC_CALL(splitAndMerge(dec,newRow,nextNode,false,newRowInformation));

        //Add other nodes in the subtree
        for (int i = 0; i < newRow->reducedMembers[nextNode].numChildren; ++i) {
            children_idx idx = newRow->reducedMembers[nextNode].firstChild + i;
            reduced_member_id child = newRow->childrenStorage[idx];
            MATREC_CALL(mergeChildrenNodes(dec,newRow,child,baseNode,newRowInformation));
        }

        //Move up one layer
        baseNode = nextNode;
        nextNode = newRow->reducedMembers[nextNode].parent;
    }

    return MATREC_OKAY;
}
static MATREC_ERROR transformComponentRowAddition(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow,
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
                MATREC_CALL(transformSingleParallel(dec,newRow,reducedMember,member,newRowInformation));
                break;
            }
            case SPQR_MEMBERTYPE_LOOP:
            case SPQR_MEMBERTYPE_SERIES: {

                newRowInformation->member = member;
                cut_arc_id cutArc = newRow->reducedMembers[reducedMember].firstCutArc;
                spqr_arc arc = newRow->cutArcs[cutArc].arc;
                newRowInformation->reversed = arcIsReversedNonRigid(dec,arc) == newRow->cutArcs[cutArc].arcReversed;
                if(type == SPQR_MEMBERTYPE_LOOP ){
                    if(getNumMemberArcs(dec,member) == 2){
                        changeLoopToSeries(dec,member);
                    }else{

                    }
                }
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


MATREC_ERROR MATRECcreateNetworkRowAddition(MATREC* env, MATRECNetworkRowAddition** pNewRow ){
    assert(env);
    MATREC_CALL(MATRECallocBlock(env,pNewRow));
    MATRECNetworkRowAddition * newRow = *pNewRow;

    newRow->remainsNetwork = true;

    newRow->reducedMembers = NULL;
    newRow->memReducedMembers = 0;
    newRow->numReducedMembers = 0;

    newRow->reducedComponents = NULL;
    newRow->memReducedComponents = 0;
    newRow->numReducedComponents = 0;

    newRow->memberInformation = NULL;
    newRow->memMemberInformation = 0;
    newRow->numMemberInformation = 0;

    newRow->cutArcs = NULL;
    newRow->memCutArcs = 0;
    newRow->numCutArcs = 0;
    newRow->firstOverallCutArc = INVALID_CUT_ARC;

    newRow->childrenStorage = NULL;
    newRow->memChildrenStorage = 0;
    newRow->numChildrenStorage = 0;

    newRow->newRowIndex = MATREC_INVALID_ROW;

    newRow->newColumnArcs = NULL;
    newRow->newColumnReversed = NULL;
    newRow->memColumnArcs = 0;
    newRow->numColumnArcs = 0;

    newRow->leafMembers = NULL;
    newRow->numLeafMembers = 0;
    newRow->memLeafMembers = 0;

    newRow->decompositionColumnArcs = NULL;
    newRow->decompositionColumnArcReversed = NULL;
    newRow->memDecompositionColumnArcs = 0;
    newRow->numDecompositionColumnArcs = 0;

    newRow->isArcCut = NULL;
    newRow->isArcCutReversed = NULL;
    newRow->memIsArcCut = 0;
    newRow->numIsArcCut = 0;

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

void MATRECfreeNetworkRowAddition(MATREC* env, MATRECNetworkRowAddition ** pNewRow){
    assert(*pNewRow);

    MATRECNetworkRowAddition * newRow = *pNewRow;
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
    MATRECfreeBlockArray(env,&newRow->isArcCut);
    MATRECfreeBlockArray(env,&newRow->isArcCutReversed);
    if(newRow->decompositionColumnArcs){
        MATRECfreeBlockArray(env,&newRow->decompositionColumnArcs);
    }
    MATRECfreeBlockArray(env,&newRow->decompositionColumnArcReversed);
    if(newRow->newColumnArcs){
        MATRECfreeBlockArray(env,&newRow->newColumnArcs);
    }
    MATRECfreeBlockArray(env,&newRow->newColumnReversed);
    if(newRow->childrenStorage){
        MATRECfreeBlockArray(env,&newRow->childrenStorage);
    }
    if(newRow->cutArcs){
        MATRECfreeBlockArray(env,&newRow->cutArcs);
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

MATREC_ERROR MATRECNetworkRowAdditionCheck(MATRECNetworkDecomposition * dec, MATRECNetworkRowAddition * newRow,
                                           const MATREC_row row, const MATREC_col * columns, const double * columnValues,
                                           size_t numColumns){
    assert(dec);
    assert(newRow);
    assert(numColumns == 0 || columns );

    newRow->remainsNetwork = true;
    cleanUpPreviousIteration(dec,newRow);

    MATREC_CALL(newRowUpdateRowInformation(dec,newRow,row,columns,columnValues,numColumns));
    MATREC_CALL(constructRowReducedDecomposition(dec,newRow));
    MATREC_CALL(createReducedDecompositionCutArcs(dec,newRow));

    MATREC_CALL(determineLeafReducedMembers(dec,newRow));
    MATREC_CALL(allocateRigidSearchMemory(dec,newRow));
    MATREC_CALL(allocateTreeSearchMemory(dec,newRow));
    //Check for each component if the cut arcs propagate through a row tree marker to a cut arc in another component
    //From the leafs inward.
    propagateComponents(dec,newRow);
    //It can happen that we are not graphic by some of the checked components.
    //In that case, further checking may lead to errors as some invariants that the code assumes will be broken.
    if(newRow->remainsNetwork){
        for (int i = 0; i < newRow->numReducedComponents; ++i) {
            determineMergeableTypes(dec,newRow,newRow->reducedComponents[i].root);
            //exit early if one is not graphic
            if(!newRow->remainsNetwork){
                break;
            }
        }
    }

    cleanUpRowMemberInformation(newRow);

    return MATREC_OKAY;
}

MATREC_ERROR MATRECNetworkRowAdditionAdd(MATRECNetworkDecomposition *dec, MATRECNetworkRowAddition *newRow){
    assert(newRow->remainsNetwork);
    if(newRow->numReducedComponents == 0){
        spqr_member newMember = SPQR_INVALID_MEMBER;
        MATREC_CALL(createStandaloneParallel(dec,newRow->newColumnArcs, newRow->newColumnReversed,
                                           newRow->numColumnArcs,newRow->newRowIndex,&newMember));
    }else if (newRow->numReducedComponents == 1){
        NewRowInformation information = emptyNewRowInformation();
        MATREC_CALL(transformComponentRowAddition(dec,newRow,&newRow->reducedComponents[0],&information));

        if(newRow->numColumnArcs == 0){
            spqr_arc rowArc = SPQR_INVALID_ARC;
            MATREC_CALL(createRowArc(dec,information.member,&rowArc,newRow->newRowIndex,information.reversed));
            if(SPQRnodeIsValid(information.head)){
                assert(SPQRnodeIsValid(information.tail));
                assert(SPQRarcIsValid(information.representative));
                //TODO: check if finds are necessary
                setArcHeadAndTail(dec,rowArc,findNode(dec,information.head),findNode(dec,information.tail));
                arcSetRepresentative(dec,rowArc,information.representative);
                arcSetReversed(dec,rowArc,information.reversed != arcIsReversedNonRigid(dec,information.representative));
            }
        }else{
            spqr_member new_row_parallel = SPQR_INVALID_MEMBER;
            MATREC_CALL(createConnectedParallel(dec,newRow->newColumnArcs,newRow->newColumnReversed,newRow->numColumnArcs,
                                              newRow->newRowIndex,&new_row_parallel));
            spqr_arc markerArc = SPQR_INVALID_ARC;
            spqr_arc ignore = SPQR_INVALID_ARC;
            MATREC_CALL(createMarkerPairWithReferences(dec,information.member,new_row_parallel,true,
                                                     information.reversed,false,
                                                     &markerArc,&ignore));
            if(SPQRnodeIsValid(information.head)){
                assert(SPQRnodeIsValid(information.tail));
                assert(SPQRarcIsValid(information.representative));
                //TODO: check if finds are necessary
                setArcHeadAndTail(dec,markerArc,findNode(dec,information.head),findNode(dec,information.tail));
                arcSetRepresentative(dec,markerArc,information.representative);
                arcSetReversed(dec,markerArc,information.reversed != arcIsReversedNonRigid(dec,information.representative));
            }
        }
        if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
            assert(getNumMemberArcs(dec,information.member) == 2 || getNumMemberArcs(dec,information.member) == 3);
            if(getNumMemberArcs(dec,information.member) == 3){
                changeLoopToSeries(dec,information.member);
            }
        }
    }else{
#ifndef NDEBUG
        int numDecComponentsBefore = numConnectedComponents(dec);
#endif
        spqr_member new_row_parallel = SPQR_INVALID_MEMBER;
        MATREC_CALL(createConnectedParallel(dec,newRow->newColumnArcs,newRow->newColumnReversed,newRow->numColumnArcs,
                                          newRow->newRowIndex,&new_row_parallel));
        for (int i = 0; i < newRow->numReducedComponents; ++i) {
            NewRowInformation information = emptyNewRowInformation();

            MATREC_CALL(transformComponentRowAddition(dec,newRow,&newRow->reducedComponents[i],&information));
            if(getMemberType(dec,information.member) == SPQR_MEMBERTYPE_LOOP){
                assert(getNumMemberArcs(dec,information.member) == 1);
                spqr_arc arc = getFirstMemberArc(dec,information.member);
                assert(newRow->isArcCut[arc]);
                moveArcToNewMember(dec, arc,information.member,new_row_parallel);
                arcSetReversed(dec,arc,newRow->isArcCutReversed[arc]);
                dec->members[information.member].type = SPQR_MEMBERTYPE_UNASSIGNED;
            }else{
                reorderComponent(dec,information.member); //Make sure the new component is the root of the local decomposition tree
                spqr_arc markerArc = SPQR_INVALID_ARC;
                spqr_arc ignore = SPQR_INVALID_ARC;
                MATREC_CALL(createMarkerPairWithReferences(dec,new_row_parallel,information.member,false,
                                                         false,information.reversed,
                                                         &ignore,&markerArc));
                if(SPQRnodeIsValid(information.head)){
                    assert(SPQRnodeIsValid(information.tail));
                    assert(SPQRarcIsValid(information.representative));
                    //TODO: check if finds are necessary
                    setArcHeadAndTail(dec,markerArc,findNode(dec,information.head),findNode(dec,information.tail));
                    arcSetRepresentative(dec,markerArc,information.representative);
                    arcSetReversed(dec,markerArc,information.reversed != arcIsReversedNonRigid(dec,information.representative));
                }
            }
        }
        decreaseNumConnectedComponents(dec,newRow->numReducedComponents-1);
        assert(numConnectedComponents(dec) == (numDecComponentsBefore - newRow->numReducedComponents + 1));
    }
    return MATREC_OKAY;
}

bool MATRECNetworkRowAdditionRemainsNetwork(const MATRECNetworkRowAddition *newRow){
    return newRow->remainsNetwork;
}
